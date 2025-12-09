#include "BlockPreconditioner.hpp"
#include "../../../Utilities/SparseMatFunctions.hpp"
#include "../../CardiovascularBCs/SumOp.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

BlockPreconditioner::BlockPreconditioner(Array<int> &offsets_,Array<FiniteElementSpace *> &fes, bool equal_order_) 
: Solver(offsets_.Last()), equal_order(equal_order_)
{
   fes.Copy(spaces);
   offsets.MakeRef(offsets_);

   Ad.SetSize(spaces[0]->GetTrueVSize());
   Ad = 0.0;
}

void BlockPreconditioner::SetOperator(const mfem::Operator &op)
{
   jacobian = (BlockOperator *) &op;

   // Setup A block solver
   SetASolver(jacobian);

   // Setup Schur block solver
   SetSchurSolver(jacobian);
   FormSchurOperator(jacobian);
}

void BlockPreconditioner::SetASolver(BlockOperator *jacobian)
{
   if (AInv == nullptr)
   {
      if (direct_solvers)
      {
         #ifdef MFEM_USE_MPI
            #ifdef MFEM_USE_MKL_CPARDISO
            CPardisoSolver *A_solver_direct = new CPardisoSolver(dynamic_cast<ParFiniteElementSpace*>(spaces[0])->GetComm());
            A_solver_direct->SetMatrixType(CPardisoSolver::MatType::REAL_NONSYMMETRIC);
            AInv = A_solver_direct;
            #else
            MFEM_ABORT("Error: Parallel direct solvers require mfem with MKL CPARDISO!");
            #endif
         #else
            #ifdef MFEM_USE_SUITESPARSE
            UMFPackSolver *A_solver_direct = new UMFPackSolver();
            A_solver_direct->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
            AInv = A_solver_direct;
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
            HypreSmoother *A_prec_jacobi = new HypreSmoother();
            A_prec_jacobi->SetType(HypreSmoother::Jacobi, 1);
            #else
            DSmoother *A_prec_jacobi = new DSmoother();
            #endif
            AInvPC = A_prec_jacobi;
         }
         else // AMG
         {
            #ifdef MFEM_USE_MPI
            HypreBoomerAMG *A_prec_amg = new HypreBoomerAMG();
            A_prec_amg->SetPrintLevel(0);
            A_prec_amg->SetSystemsOptions(spaces[0]->GetMesh()->Dimension(),spaces[0]->GetOrdering() == Ordering::byNODES ? true : false); // RB 2-7-24
            
            // Additonal AMG options
            int dim = spaces[0]->GetMesh()->Dimension();
            // Cycle type
            A_prec_amg->SetCycleType(1); // 1 = V-cycle, 2 = W-cycle
            // Coarsening method
            A_prec_amg->SetCoarsening(10);   // 10 = HMIS, 8 = PMIS, 6 = Falgout, 0 = CLJP)
            // Interp method
            A_prec_amg->SetInterpolation(6);    // 6 = extended+i, 0 = classical)
            // Smoother/Relaxation type
            A_prec_amg->SetRelaxType(3); // hybrid gauss-seidel smoother, 18 = l1 jacobi MFEM default
            // // Truncation factor for interp
            // HYPRE_BoomerAMGSetTruncFactor(A_prec_amg,0.3);
            // Strong threshold for 3D is 0.7, for 2D is 0.25
            A_prec_amg->SetStrengthThresh(dim == 3 ? 0.7 : 0.25);
            // // Max number of elements per row for interp
            // HYPRE_BoomerAMGSetPMaxElmts(A_prec_amg, 5);
            // Set aggressive coarsening levels (default 1)
            A_prec_amg->SetAggressiveCoarsening(2);
            
            // A_prec_amg->SetRelaxType(3); // non symm smoother for non symmetric A matrix // RB 11/15/23
            
            AInvPC = A_prec_amg;
            #else
            MFEM_ABORT("Error: AMG requires mfem with MPI!");
            #endif
         }

         GMRESSolver *A_solver_iter = nullptr;
         #ifdef MFEM_USE_MPI
         if (!svsolver && cv_bcs) A_solver_iter = new ModifiedGMRESSolver(dynamic_cast<ParFiniteElementSpace *>(spaces[0])->GetComm());
         else A_solver_iter = new GMRESSolver(dynamic_cast<ParFiniteElementSpace *>(spaces[0])->GetComm());
         #else
         if (!svsolver) A_solver_iter = new ModifiedGMRESSolver();
         else A_solver_iter = new GMRESSolver();
         #endif
         A_solver_iter->SetRelTol(a_rtol);
         A_solver_iter->SetAbsTol(a_atol);
         A_solver_iter->SetMaxIter(a_maxIt);
         // A_solver_iter->SetKDim(a_maxIt);
         A_solver_iter->SetPrintLevel(-1);
         if (!svsolver && cv_bcs) A_solver_iter->SetPreconditioner(*AInvPC);  
         else A_solver_iter->SetPreconditioner(cv_bcs ? *cvbc_PC : *AInvPC);     
         A_solver_iter->iterative_mode = false;

         AInv = A_solver_iter;
      }
   }
   // At each Newton cycle, compute the A block preconditioner by
   // updating the iterative solver which, in turn, updates its preconditioner
   if (cv_bcs && !svsolver)
   {
      // form sparse approx of A = diag(Abc) + A3d
      Vector &Abc_diag = cvbc_PC->GetAbcDiag();
      
      const SumOp *sum_op = dynamic_cast<const SumOp *>(&jacobian->GetBlock(0,0));
      #ifdef MFEM_USE_MPI
      HypreParMatrix *A3d = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(sum_op->GetA()));
      A_cvbc = new HypreParMatrix(*A3d); // A_cvbc is owned by this class
      HypreParMatrixAddToDiag(A_cvbc,Abc_diag); // add Abc_diag to diag of A3d
      #else
      SparseMatrix *A3d = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(sum_op->GetA()));
      A_cvbc = new SparseMatrix(*A3d); // A_cvbc is owned by this class
      SparseMatAddToDiag(*A_cvbc,Abc_diag); // add Abc_diag to diag of A3d
      #endif

      dynamic_cast<ModifiedGMRESSolver *>(AInv)->SetOperatorGMRES(jacobian->GetBlock(0,0));
      dynamic_cast<ModifiedGMRESSolver *>(AInv)->SetOperatorPrec(*A_cvbc);
   }
   else
   {
      AInv->SetOperator(jacobian->GetBlock(0,0));
   }
}

void BlockPreconditioner::SetSchurSolver(BlockOperator *jacobian)
{
   // setup Schur complement operator and its inverse solver
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
            int dim = spaces[0]->GetMesh()->Dimension();
            // Cycle type
            S_prec_amg->SetCycleType(1); // 1 = V-cycle, 2 = W-cycle
            // Coarsening method
            S_prec_amg->SetCoarsening(10);   // 10 = HMIS, 8 = PMIS, 6 = Falgout, 0 = CLJP)
            // Interp method
            S_prec_amg->SetInterpolation(6);    // 6 = extended+i, 0 = classical)
            // Smoother/Relaxation type
            S_prec_amg->SetRelaxType(3); // hybrid gauss-seidel smoother, 18 = l1 jacobi MFEM default
            // // Truncation factor for interp
            // HYPRE_BoomerAMGESetTruncFactor(S_prec_amg,0.3);
            // Strong threshold for 3D is 0.7, for 2D is 0.25
            S_prec_amg->SetStrengthThresh(dim == 3 ? 0.7 : 0.25);
            // // Max number of elements per row for interp
            // HYPRE_BoomerAMGSetPMaxElmts(S_prec_amg, 5);
            // Set aggressive coarsening levels (default 1)
            S_prec_amg->SetAggressiveCoarsening(2);


            SInvPC = S_prec_amg;
            #else
            MFEM_ABORT("Error: AMG requires mfem with MPI!");
            #endif
         }

         #ifdef MFEM_USE_MPI
         GMRESSolver *S_solver_iter = new GMRESSolver(dynamic_cast<ParFiniteElementSpace *>(spaces[0])->GetComm());
         #else
         GMRESSolver *S_solver_iter = new GMRESSolver();
         #endif
         S_solver_iter->SetRelTol(s_rtol);
         S_solver_iter->SetAbsTol(s_atol);
         S_solver_iter->SetMaxIter(s_maxIt);
         // S_solver_iter->SetKDim(s_maxIt);
         S_solver_iter->SetPreconditioner(*SInvPC);
         S_solver_iter->iterative_mode = false;

         SInv = S_solver_iter;
      }      
   }
}

void BlockPreconditioner::MultTranspose(const Vector & x,Vector & y) const
{
   MFEM_ABORT("Error: To Do");
}

BlockPreconditioner::~BlockPreconditioner()
{
   delete S;
   delete tempS;
   delete SInv;
   delete SInvPC;
   delete AInv;
   delete AInvPC;
   delete A_cvbc;

   // set to nullptr
   jacobian = nullptr;
   tempA = nullptr;
   for (int i = 0; i < spaces.Size(); i++)
   {
      spaces[i] = nullptr;
   }
}