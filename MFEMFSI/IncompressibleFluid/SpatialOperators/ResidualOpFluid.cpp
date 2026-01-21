#include "ResidualOpFluid.hpp"
#include <general/forall.hpp>
#include <iomanip>

using namespace mfem;
using namespace std;
using namespace ifs_solver;
using namespace immersed_solid;

void ResidualOpFluid::Mult(const Vector &k, Vector &y) const
{   
   MFEM_ABORT("Not implemented");
}

void ResidualOpFluid::Mult_PC(const Vector &y, const Vector &ydot, const Vector &ft, const Vector &fb, double am, Vector &res, bool line_search) const
{
   // Extract the blocks from the input and output vectors
   rblock.Update(res.GetData(),block_trueOffsets); // [rv rp] // TO DO check if this allocates mem each time!
   Vector &rv = rblock.GetBlock(0);
   Vector &rp = rblock.GetBlock(1);

   yblock.Update(y.GetData(),block_trueOffsets); // [v p]
   Vector &v = yblock.GetBlock(0);
   Vector &p = yblock.GetBlock(1);

   ydotblock.Update(ydot.GetData(),block_trueOffsets); // [vdot pdot]
   Vector &vdot = ydotblock.GetBlock(0);
   Vector &pdot = ydotblock.GetBlock(1);
   
   if (total_integ)
   {
      sw_fluidres->Start();
      // Let integrator (FluidNF) know about current vdot
      curr_vdot->SetFromTrueDofs(vdot);

      // fluid integrals (block nonlinearform)
      FluidNF->Mult(y,res);

      // external work
      // // body forces (handled in above nonlinear form)
      // rv.Add(-1.0,fb);

      // surface tractions
      rv.Add(-1.0,ft);

      // CV bcs
      if (cv_bcs)
      {
         cvbc_op->AddMult(v,rv); // overloaded to be efficient
      }
      sw_fluidres->Stop();

      if (IFEM)
      {
         im_solid->ComputeSolidResidualAndSpread2Fluid(*f_fsi);
         sw_fluidres->Start();       
         rv.Add(1.0,*f_fsi);
      }
   }
   else
   {
      MFEM_ABORT("TO DO");
   }

   ////////////////////////////////////
   // Zero out essential dofs in RHS vector
   // I implement dirichlet bcs in the following way:
   // Set essential dofs to n+1 in the predictors
   // Solve for deltas to essential dofs that are zero

   // velocity
   const int N = vel_ess_tdofs->Size();
   const auto idx = vel_ess_tdofs->Read();
   auto rv_RW = rv.ReadWrite();
   MFEM_FORALL(i, N, rv_RW[idx[i]] = 0.0;);

   // pressure
   const int Np = pres_ess_tdofs->Size();
   const auto idxp = pres_ess_tdofs->Read();
   auto rp_RW = rp.ReadWrite();
   MFEM_FORALL(j, Np, rp_RW[idxp[j]] = 0.0;);
   sw_fluidres->Stop();
}

Operator &ResidualOpFluid::GetGradient(const Vector &k) const
{
   MFEM_ABORT("Not implemented");
   return *gradOp;
}

Operator &ResidualOpFluid::GetGradient_PC(const Vector &y, const Vector &ydot, double dt, double am, double af) const
{
   if (direct_solver || monolithic_prec)
   {
      #ifndef MFEM_USE_MPI
      if(gradMat) delete gradMat;
      gradMat = new BlockMatrix(block_trueOffsets);
      #endif
   }
   else
   {
      if(gradOp) delete gradOp;
      gradOp = new BlockOperator(block_trueOffsets);
   }
   // use exisiting data in yblock!

   // yblock.Update(y.GetData(),block_trueOffsets); // [v p]
   // Vector &v = yblock.GetBlock(0);
   // Vector &p = yblock.GetBlock(1);

   // ydotblock.Update(ydot.GetData(),block_trueOffsets); // [vdot pdot]
   // Vector &vdot = ydotblock.GetBlock(0);
   // Vector &pdot = ydotblock.GetBlock(1);

   // // get current vdot if gen alpha and stabilized method
   // curr_vdot->SetFromTrueDofs(vdot);

   if (total_integ)
   {
      // Compute the tangent matrix
      sw_fluidtan->Start();
      ComputeGradient(y);

      // Modify for essential dofs, CVBCs
      // (0,0), expects Agrad to hold the (0,0) block
      if (!cv_bcs)
      {
         #ifdef MFEM_USE_MPI
         // const_cast<HypreParMatrix *>(Agrad)->EliminateBC(*vel_ess_tdofs,DIAG_ONE);
         Agrad->EliminateBC(*vel_ess_tdofs,DIAG_ONE);
         #else
         // const_cast<SparseMatrix *>(Agrad)->EliminateBC(*vel_ess_tdofs,DIAG_ONE);
         Agrad->EliminateBC(*vel_ess_tdofs,DIAG_ONE);
         #endif
         if (direct_solver || monolithic_prec)
         {
            #ifdef MFEM_USE_MPI
            gradMat(0,0) = Agrad;
            #else
            gradMat->SetBlock(0,0,Agrad);
            #endif
         } 
         else gradOp->SetBlock(0,0,Agrad);
      }  
      else
      {
         // CV BCs
         // get diag of A (without CV bcs)
         #ifdef MFEM_USE_MPI
         // const_cast<HypreParMatrix *>(Agrad)->GetDiag(*A3d_diag);
         Agrad->GetDiag(*A3d_diag);
         #else
         // const_cast<SparseMatrix *>(Agrad)->GetDiag(*A3d_diag);
         Agrad->GetDiag(*A3d_diag);
         #endif

         // get diag of Abc
         cvbc_op->GetDiag(*Abc_diag,af); // Atotal = Abc
         
         // Atotal = Abc + A3d
         const int Nsize = A3d_diag->Size();
         const bool use_dev = A3d_diag->UseDevice();
         auto yvec = Atotal_diag->ReadWrite(use_dev);
         auto xvec = A3d_diag->Read(use_dev);
         auto xvec2 = Abc_diag->Read(use_dev);
         mfem::forall_switch(use_dev, Nsize, [=] MFEM_HOST_DEVICE (int i) { yvec[i] = xvec[i] + xvec2[i]; });

         if (sum_grad_op) delete sum_grad_op; // Do I need to delete this?
         sum_grad_op = new SumOp(vfes->GetTrueVSize(),*vel_ess_tdofs);
         sum_grad_op->SetOps(*Agrad,*cvbc_op);
         sum_grad_op->SetCV_timefactor(af);

         if (direct_solver || monolithic_prec) 
         {
            MFEM_ABORT("Direct solver not supported with CV BCs!");
         }
         else 
         {
            gradOp->SetBlock(0,0,sum_grad_op); // Note that essential dofs are applied in the action of SumOp Class
         }
      }      

      // (0,1)
      #ifdef MFEM_USE_MPI
      HypreParMatrix *temp_e = G_stab->EliminateCols(*pres_ess_tdofs);
      G_stab->EliminateRows(*vel_ess_tdofs);
      delete temp_e;
      #else
      Array<int> ess_trial_tdof_marker;
      FiniteElementSpace::ListToMarker(*pres_ess_tdofs, pfes->GetTrueVSize(),
                                       ess_trial_tdof_marker);
      SparseMatrix *mat_e = new SparseMatrix(const_cast<SparseMatrix *>(G_stab)->Height(), const_cast<SparseMatrix *>(G_stab)->Width());
      const_cast<SparseMatrix *>(G_stab)->EliminateCols(ess_trial_tdof_marker, *mat_e);
      for (int i=0; i<vel_ess_tdofs->Size(); ++i)
      {
         const_cast<SparseMatrix *>(G_stab)->EliminateRow((*vel_ess_tdofs)[i]);
      }
      delete mat_e;   
      #endif
      if (direct_solver || monolithic_prec) 
      {
         #ifdef MFEM_USE_MPI
         gradMat(0,1) = G_stab;
         #else
         gradMat->SetBlock(0,1,G_stab);
         #endif
      }
      else gradOp->SetBlock(0,1,G_stab);

      
      // (1,0)
      #ifdef MFEM_USE_MPI
      HypreParMatrix *temp_e_2 = D_stab->EliminateCols(*vel_ess_tdofs);
      D_stab->EliminateRows(*pres_ess_tdofs);
      delete temp_e_2;
      #else
      Array<int> ess_trial_tdof_marker_v;
      FiniteElementSpace::ListToMarker(*vel_ess_tdofs, vfes->GetTrueVSize(),
                                       ess_trial_tdof_marker_v);

      SparseMatrix *mat_e_2 = new SparseMatrix(const_cast<SparseMatrix *>(D_stab)->Height(), const_cast<SparseMatrix *>(D_stab)->Width());
      const_cast<SparseMatrix *>(D_stab)->EliminateCols(ess_trial_tdof_marker_v, *mat_e_2);
      for (int i=0; i<pres_ess_tdofs->Size(); ++i)
      {
         const_cast<SparseMatrix *>(D_stab)->EliminateRow((*pres_ess_tdofs)[i]);
      }
      delete mat_e_2;
      #endif
      if (direct_solver || monolithic_prec) 
      {
         #ifdef MFEM_USE_MPI
         gradMat(1,0) = D_stab;
         #else
         gradMat->SetBlock(1,0,D_stab);
         #endif
      }
      else gradOp->SetBlock(1,0,D_stab);

      if (equal_order)
      {
         #ifdef MFEM_USE_MPI
         C_stab->EliminateBC(*pres_ess_tdofs,DIAG_ONE);
         #else
         const_cast<SparseMatrix *>(C_stab)->EliminateBC(*pres_ess_tdofs,DIAG_ONE);
         #endif
         if (direct_solver || monolithic_prec) 
         {
            #ifdef MFEM_USE_MPI
            gradMat(1,1) = C_stab;
            #else
            gradMat->SetBlock(1,1,C_stab);
            #endif
         }
         else
         {
            gradOp->SetBlock(1,1,C_stab);
         } 
      }      
   }
   else
   {
      MFEM_ABORT("TO DO");
   }

   if (direct_solver || monolithic_prec)
   {
      if (totalGradMat) delete totalGradMat;  

      #ifdef MFEM_USE_MPI
      totalGradMat = HypreParMatrixFromBlocks(gradMat);
      #else
      totalGradMat = dynamic_cast<BlockMatrix *>(gradMat)->CreateMonolithic();
      #endif

      sw_fluidtan->Stop();
      return *totalGradMat;
   }
   else
   {
      sw_fluidtan->Stop();
      return *gradOp;
   }
}

#ifdef MFEM_USE_MPI
void ResidualOpFluid::ComputeGradient(const Vector &y) const
{
   BlockOperator *grad_stab = dynamic_cast<BlockOperator *>(&FluidNF->GetGradient(y));
   A_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(0,0)));
   G_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(0,1)));
   D_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(1,0)));
   if (equal_order) C_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(1,1)));

   if (IFEM)
   {
      sw_fluidtan->Stop();
      if (Agrad) delete Agrad;
      // get C K_solid Ct and add to (0,0) block
      HypreParMatrix *CKCt = dynamic_cast<HypreParMatrix *>(&im_solid->FormLHS());
      sw_fluidtan->Start();
      Agrad = ParAdd(A_stab,CKCt);
   }
   else
   {
      Agrad = A_stab;
   }
}
#else
void ResidualOpFluid::ComputeGradient(const Vector &y) const
{
   #ifdef WIN32
   BlockOperator *grad_stab = dynamic_cast<BlockOperator *>(&FluidNF->GetGradient(y));
   A_stab = dynamic_cast<SparseMatrix *>(&grad_stab->GetBlock(0,0));
   G_stab = dynamic_cast<SparseMatrix *>(&grad_stab->GetBlock(0,1));
   D_stab = dynamic_cast<SparseMatrix *>(&grad_stab->GetBlock(1,0));  
   if (equal_order) C_stab = dynamic_cast<SparseMatrix *>(&grad_stab->GetBlock(1,1));
   #else
   BlockOperator *grad_stab = dynamic_cast<BlockOperator *>(&FluidNF->GetGradient(y));
   A_stab = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&grad_stab->GetBlock(0,0)));
   G_stab = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&grad_stab->GetBlock(0,1)));
   D_stab = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&grad_stab->GetBlock(1,0)));
   if (equal_order) C_stab = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&grad_stab->GetBlock(1,1)));
   #endif

   // (0,0)
   if (IFEM)
   {
      sw_fluidtan->Stop();
      if (Agrad) delete Agrad;
      // get C^t K_solid C and add to grad
      SparseMatrix *CtKC = dynamic_cast<SparseMatrix *>(&im_solid->FormLHS());
      sw_fluidtan->Start();
      Agrad = Add(*A_stab,*CtKC);
   }
   else
   {
      Agrad = A_stab;
   }
}
#endif

ResidualOpFluid::~ResidualOpFluid()
{
   delete gradOp;
   delete sum_grad_op;

   if (IFEM) delete Agrad;

   delete totalGradMat;
   #ifdef MFEM_USE_MPI
   gradMat = nullptr;
   #else   
   delete gradMat;
   #endif
}