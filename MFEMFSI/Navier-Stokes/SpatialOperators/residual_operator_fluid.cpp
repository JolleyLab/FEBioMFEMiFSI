#include "residual_operator_fluid.hpp"
#include <general/forall.hpp>
#include <iomanip>

using namespace mfem;
using namespace std;
using namespace ins_solver;
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

   // double norm = mfem::InnerProduct(MPI_COMM_WORLD,rblock,rblock);
   // // double norm = mfem::InnerProduct(rblock,rblock);
   // mfem::out << "Residual norm = " << std::setprecision(16) << sqrt(norm) << std::endl;
   // // MFEM_ABORT("HERE");

   // // print out residual vector to file
   // std::ofstream outfile("r_vec.txt", std::ios::out);
   // int ii,jj;
   // int precision_out = 6;
   // int width_out = precision_out + 6;
   // for (int i = 0; i < rblock.Size(); ++i) 
   // {
   //    if (abs(rblock(i)) >= 1e-12) outfile << std::setw(width_out) << std::setprecision(precision_out) << rblock(i);
   //    else outfile << std::setw(width_out) << std::setprecision(precision_out) << 0.0;
   //    outfile << std::endl;
   // }
   // outfile << std::endl;
   
   // // Close the file
   // outfile.close();

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
      // #ifndef MFEM_USE_MPI
      // #ifdef WIN32
      // BlockOperator *grad_stab = dynamic_cast<BlockOperator *>(&FluidNF->GetGradient(y));
      // SparseMatrix *A_stab = dynamic_cast<SparseMatrix *>(&grad_stab->GetBlock(0,0));
      // SparseMatrix *G_stab = dynamic_cast<SparseMatrix *>(&grad_stab->GetBlock(0,1));
      // SparseMatrix *D_stab = dynamic_cast<SparseMatrix *>(&grad_stab->GetBlock(1,0));  
      // #else
      // BlockOperator *grad_stab = dynamic_cast<BlockOperator *>(&FluidNF->GetGradient(y));
      // SparseMatrix *A_stab = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&grad_stab->GetBlock(0,0)));
      // SparseMatrix *G_stab = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&grad_stab->GetBlock(0,1)));
      // SparseMatrix *D_stab = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&grad_stab->GetBlock(1,0)));
      // #endif
      // #else
      // BlockOperator *grad_stab = dynamic_cast<BlockOperator *>(&FluidNF->GetGradient(y));
      // HypreParMatrix *A_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(0,0)));
      // HypreParMatrix *G_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(0,1)));
      // HypreParMatrix *D_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(1,0)));
      // #endif

      // // Modify for essential dofs, CVBCs, immersed solids
      // // (0,0)
      // if (IFEM)
      // {
      //    #ifdef MFEM_USE_MPI
      //    MFEM_ABORT("TODO");
      //    // need to return hypreparmatrix from immersed solid even if serial
      //    if (Agrad) delete Agrad;
      //    // get C^t K_solid C and add to grad
      //    HypreParMatrix *CtKC = dynamic_cast<HypreParMatrix *>(&im_solid->FormLHS());
      //    Agrad = ParAdd(A_stab,CtKC);
      //    #else
      //    if (Agrad) delete Agrad;
      //    // get C^t K_solid C and add to grad
      //    SparseMatrix *CtKC = dynamic_cast<SparseMatrix *>(&im_solid->FormLHS());
      //    Agrad = Add(*A_stab,*CtKC);
      //    #endif
      // }
      // else
      // {
      //    Agrad = A_stab;
      // }


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
      // HypreParMatrix *temp_e = const_cast<HypreParMatrix *>(G_stab)->EliminateCols(*pres_ess_tdofs);
      // const_cast<HypreParMatrix *>(G_stab)->EliminateRows(*vel_ess_tdofs);
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
      // HypreParMatrix *temp_e_2 = const_cast<HypreParMatrix *>(D_stab)->EliminateCols(*vel_ess_tdofs);
      // const_cast<HypreParMatrix *>(D_stab)->EliminateRows(*pres_ess_tdofs);
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
         // HypreParMatrix *C_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(1,1)));
         // const_cast<HypreParMatrix *>(C_stab)->EliminateBC(*pres_ess_tdofs,DIAG_ONE);
         C_stab->EliminateBC(*pres_ess_tdofs,DIAG_ONE);
         #else
         // #ifdef WIN32
         // SparseMatrix *C_stab = dynamic_cast<SparseMatrix *>(&grad_stab->GetBlock(1,1));   
         // #else
         // SparseMatrix *C_stab = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&grad_stab->GetBlock(1,1))); 
         // #endif    
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

   // // // Print gradOp
   // // std::ofstream A00mat_file("A00mat.txt");
   // std::ofstream A00mat_file("A00mat_febio.txt");
   // SparseMatrix* A00 = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&gradMat->GetBlock(0,0)));
   // if (A00mat_file.is_open())
   // {
   //    A00->PrintMatlab(A00mat_file);
   //    A00mat_file.close();
   // }

   // // Print gradOp block (0,1)
   // // std::ofstream A01mat_file("A01mat.txt");
   // std::ofstream A01mat_file("A01mat_v2.txt");
   // SparseMatrix* A01 = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&gradMat->GetBlock(0,1)));
   // if (A01mat_file.is_open())
   // {
   //    A01->PrintMatlab(A01mat_file);
   //    A01mat_file.close();
   // }

   // // Print gradOp block (1,0)
   // // std::ofstream A10mat_file("A10mat.txt");
   // std::ofstream A10mat_file("A10mat_v2.txt");
   // SparseMatrix* A10 = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&gradMat->GetBlock(1,0)));
   // if (A10mat_file.is_open())
   // {
   //    A10->PrintMatlab(A10mat_file);
   //    A10mat_file.close();
   // }
   
   // // Print gradOp block (1,1)
   // // std::ofstream A11mat_file("A11mat.txt");
   // std::ofstream A11mat_file("A11mat_v2.txt");
   // SparseMatrix* A11 = const_cast<SparseMatrix *>(dynamic_cast<const SparseMatrix *>(&gradMat->GetBlock(1,1)));
   // if (A11mat_file.is_open())
   // {
   //    A11->PrintMatlab(A11mat_file);
   //    A11mat_file.close();
   // }

   if (direct_solver || monolithic_prec)
   {
      // #ifdef MFEM_USE_MPI
      // if (totalGradMat) delete totalGradMat;
      // totalGradMat = HypreParMatrixFromBlocks(gradMat);

      // sw_fluidtan->Stop();
      // return *totalGradMat;
      // #else
      // if (totalGradMat) delete totalGradMat;   
      // totalGradMat = dynamic_cast<BlockMatrix *>(gradMat)->CreateMonolithic();

      // sw_fluidtan->Stop();
      // return *totalGradMat;
      // #endif

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

      // // string file_prefix = "my_matA_mfem";
      // string file_prefix = "my_matA_febio";
      // CKCt->Print(file_prefix.c_str()); 
      // MFEM_ABORT("HERE");
   }
   else
   {
      Agrad = A_stab;
   }

   // string file_prefix = "my_matA_mfem";
   // // string file_prefix = "my_matA_febio";
   // Agrad->Print(file_prefix.c_str()); 
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