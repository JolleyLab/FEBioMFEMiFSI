#include "SCR_Prec.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

#ifdef MFEM_USE_MPI
SCR_Prec::SCR_Prec(Array<int> &offsets_,Array<FiniteElementSpace *> &fes, bool equal_order_) 
: BlockPreconditioner(offsets_,fes,equal_order_)
{
   bhatp.SetSize(offsets[2]-offsets[1]);
   bhatp = 0.0;
   yhatu.SetSize(offsets[1]-offsets[0]);
   yhatu = 0.0;

   // resize temp vectors
   tmp1.SetSize(offsets[2]-offsets[1]);
   tmp1 = 0.0;
   tmp2.SetSize(offsets[1]-offsets[0]);
   tmp2 = 0.0;
   tmp3.SetSize(offsets[1]-offsets[0]);
   tmp3 = 0.0;
   
   schur_comp = new SchurOp(offsets[2]-offsets[1],offsets[1]-offsets[0],equal_order_);
}

void SCR_Prec::SetASolver(mfem::BlockOperator *jacobian)
{
   BlockPreconditioner::SetASolver(jacobian);

   if (dynamic_cast<ModifiedGMRESSolver*>(AInv))
   {
      if (A_inner_Inv == nullptr)
      {
         auto A_solver_type = dynamic_cast<ModifiedGMRESSolver*>(AInv);
         A_inner_Inv = new ModifiedGMRESSolver(*A_solver_type); // deep copy 

         // reset tolerances for inner AInv solve of matrix-free schur complement operator
         dynamic_cast<IterativeSolver *>(A_inner_Inv)->SetRelTol(a_inner_rtol);
         dynamic_cast<IterativeSolver *>(A_inner_Inv)->SetAbsTol(a_inner_atol);
      } 
      dynamic_cast<ModifiedGMRESSolver *>(A_inner_Inv)->SetOperatorGMRES(jacobian->GetBlock(0,0));
      dynamic_cast<ModifiedGMRESSolver *>(A_inner_Inv)->SetOperatorPrec(*A_cvbc);
   }
   else
   {
      if(A_inner_Inv == nullptr)
      {
         auto A_solver_type = dynamic_cast<GMRESSolver*>(AInv);
         A_inner_Inv = new GMRESSolver(*A_solver_type); // deep copy

         // reset tolerances for inner AInv solve of matrix-free schur complement operator
         dynamic_cast<IterativeSolver *>(A_inner_Inv)->SetRelTol(a_inner_rtol);
         dynamic_cast<IterativeSolver *>(A_inner_Inv)->SetAbsTol(a_inner_atol);
      }
      A_inner_Inv->SetOperator(jacobian->GetBlock(0,0));
   }
}

void SCR_Prec::SetSchurSolver(mfem::BlockOperator *jacobian)
{
   if (SInv == nullptr)
   {
      if (direct_solvers)
      {  
         #ifdef MFEM_USE_MPI
            #ifdef MFEM_USE_MKL_CPARDISO
            CPardisoSolver *S_solver_direct = new CPardisoSolver(dynamic_cast<ParFiniteElementSpace*>(spaces[1])->GetComm());
            S_solver_direct->SetMatrixType(CPardisoSolver::MatType::REAL_NONSYMMETRIC);
            SInv = S_solver_direct;
            #else
            MFEM_ABORT("Error: Parallel direct solvers require mfem with MKL CPARDISO!");
            #endif
         #else
            #ifdef MFEM_USE_SUITESPARSE
            UMFPackSolver *S_solver_direct = new UMFPackSolver();
            S_solver_direct->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
            SInv = S_solver_direct;
            #else
            MFEM_ABORT("Error: Serial direct solvers require mfem with SuiteSparse!");
            #endif
         #endif
      }
      else
      {
         if (jacobi_preconditioning)
         {
            #ifdef MFEM_USE_MPI
            HypreSmoother *S_prec_jacobi = new HypreSmoother();
            S_prec_jacobi->SetType(HypreSmoother::Jacobi, 1);
            #else
            DSmoother *S_prec_jacobi = new DSmoother();
            #endif
            SInvPC = S_prec_jacobi;
         }
         else // AMG
         {
            #ifdef MFEM_USE_MPI
            HypreBoomerAMG *S_prec_amg = new HypreBoomerAMG();
            S_prec_amg->SetPrintLevel(0);
            S_prec_amg->SetSystemsOptions(spaces[0]->GetMesh()->Dimension(),spaces[0]->GetOrdering() == Ordering::byNODES ? true : false);
            
            // Additonal AMG options
            // Cycle type
            S_prec_amg->SetCycleType(amg_set.cycle_type); // 1 = V-cycle, 2 = W-cycle
            // Coarsening method
            S_prec_amg->SetCoarsening(amg_set.coarsening);   // 10 = HMIS, 8 = PMIS, 6 = Falgout, 0 = CLJP)
            // Interp method
            S_prec_amg->SetInterpolation(amg_set.interp);    // 6 = extended+i, 0 = classical)
            // Smoother/Relaxation type
            S_prec_amg->SetRelaxType(amg_set.relax_type); // hybrid gauss-seidel smoother, 18 = l1 jacobi MFEM default
            // // Truncation factor for interp
            // HYPRE_BoomerAMGESetTruncFactor(S_prec_amg,0.3);
            // Strong threshold for 3D is 0.7, for 2D is 0.25
            S_prec_amg->SetStrengthThresh(amg_set.strength_thresh);
            // // Max number of elements per row for interp
            // HYPRE_BoomerAMGSetPMaxElmts(S_prec_amg, 5);
            // Set aggressive coarsening levels (default 1)
            S_prec_amg->SetAggressiveCoarsening(amg_set.aggressive_coarsening);

            SInvPC = S_prec_amg;
            #else
            MFEM_ABORT("Error: AMG requires mfem with MPI!");
            #endif
         }

         ModifiedGMRESSolver *S_solver_iter = new ModifiedGMRESSolver(dynamic_cast<ParFiniteElementSpace *>(spaces[1])->GetComm());
         S_solver_iter->SetRelTol(s_rtol);
         S_solver_iter->SetAbsTol(s_atol);
         S_solver_iter->SetMaxIter(s_maxIt);
         // S_solver_iter->SetKDim(s_maxIt);
         S_solver_iter->SetPrintLevel(-1);
         S_solver_iter->SetPreconditioner(*SInvPC);
         S_solver_iter->iterative_mode = false;

         // use matrix free Schur Complement in GMRES iterations
         // but construct amg preconditioner using sparse approx of S below
         // schur_comp->SetAinvSolver(AInv);
         schur_comp->SetAinvSolver(A_inner_Inv);
         S_solver_iter->SetOperatorGMRES(*schur_comp); 

         SInv = S_solver_iter;
      }      
   }
}

void SCR_Prec::FormSchurOperator(BlockOperator *jacobian)
{
   // form Sparse approx of S = C - D(diagA)^-1 * G, formed each Newton iteration
   if (cv_bcs)
   {
      cvbc_PC->GetTotalDiag(Ad);
   }
   else
   {
      tempA = (HypreParMatrix *) &jacobian->GetBlock(0,0);
      tempA->GetDiag(Ad);
   }

   if(tempS) delete tempS;
   // any coeffs on G and D blocks, all other blocks will have correct factors already
   Gcoeff = jacobian->GetBlockCoef(0,1);
   Dcoeff = jacobian->GetBlockCoef(1,0);

   tempG = (HypreParMatrix *) &jacobian->GetBlock(0,1);
   HypreParMatrixInvScaleRows(tempG,Ad,Gcoeff*Dcoeff*-1.0); // G = G_0 * 1/diag(A) * (af * af * -1.0)

   tempD = (HypreParMatrix *) &jacobian->GetBlock(1,0);
   tempS = ParMult(tempD,tempG,true);

   HypreParMatrixScaleRows(tempG,Ad,1.0/(Gcoeff*Dcoeff*-1.0)); // G = [G_0 * 1/diag(A) * (af * af * -1.0)] * [diag(A) / (af * af * -1.0)] = G_0

   if (equal_order)
   {
      if (S) delete S;
      
      tempC = (HypreParMatrix *) &jacobian->GetBlock(1,1);
      S = ParAdd(tempS,tempC);
      // SInv->SetOperator(*S);  
      static_cast<ModifiedGMRESSolver *>(SInv)->SetOperatorPrec(*S);   
   }
   else
   {
      // SInv->SetOperator(*tempS); 
      static_cast<ModifiedGMRESSolver *>(SInv)->SetOperatorPrec(*tempS); 
   }

   // static_cast<HypreParMatrix *>(tempS)->EliminateZeroRows(); // this is a little bit faster than above!
   if (equal_order)
   {
      S->EliminateZeroRows();
   }
   else
   {
      tempS->EliminateZeroRows();
   }

   // form matrix free schur complement operator
   schur_comp->SetOperators(*tempD,Dcoeff,*tempG,Gcoeff,*tempC);
}

// Operator application
void SCR_Prec::Mult(const Vector & x,Vector & y) const
{
   MFEM_ASSERT(x.Size() == width, "incorrect input Vector size");
   MFEM_ASSERT(y.Size() == height, "incorrect output Vector size");

   yblock.Update(y.GetData(),offsets); // [u p]
   xblock.Update(x.GetData(),offsets); // [bu bp]

   // Step 1
   AInv->Mult(xblock.GetBlock(0),yhatu);
   tempD->Mult(Dcoeff,yhatu,0.0,tmp1);
   add(1.0,xblock.GetBlock(1),-1.0,tmp1,bhatp);

   // Step 2
   SInv->Mult(bhatp,yblock.GetBlock(1));

   // Step 3
   tempG->Mult(Gcoeff,yblock.GetBlock(1),0.0,tmp2);
   AInv->Mult(tmp2,tmp3);
   add(1.0,yhatu,-1.0,tmp3,yblock.GetBlock(0));
}

SCR_Prec::~SCR_Prec()
{
   delete schur_comp;
   delete A_inner_Inv;
}
#endif