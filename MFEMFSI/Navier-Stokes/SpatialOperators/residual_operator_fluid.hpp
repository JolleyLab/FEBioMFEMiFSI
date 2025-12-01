#pragma once

#include "mfem.hpp"
#include "../../Utilities/hypre_func_extensions.hpp"
#include "../../Utilities/extended_operators.hpp"
#include "../CardiovascularBCs/CV_Op.hpp"
#include "../../ImmersedFSI/immersed_solid.hpp"

using namespace immersed_solid;

namespace ins_solver
{

// class for nonlinear residual operator for newton solver
// implements computation of the residual and jacobian
class ResidualOpFluid : public mfem::Operator
{   
private:
   // temps
   mutable mfem::BlockVector rblock, yblock, ydotblock;
   mfem::Array<int> block_trueOffsets;
   mutable mfem::Vector temp; // size = u + p size
   mutable mfem::BlockOperator *gradOp = nullptr; // owned
   mutable mfem::Vector f_fsi_loc;
   #ifdef MFEM_USE_MPI
   mutable mfem::HypreParMatrix *A_stab = nullptr; // not owned
   mutable mfem::HypreParMatrix *G_stab = nullptr; // not owned
   mutable mfem::HypreParMatrix *D_stab = nullptr; // not owned
   mutable mfem::HypreParMatrix *C_stab = nullptr; // not owned
   #else
   mutable mfem::SparseMatrix *A_stab = nullptr; // not owned
   mutable mfem::SparseMatrix *G_stab = nullptr; // not owned
   mutable mfem::SparseMatrix *D_stab = nullptr; // not owned
   mutable mfem::SparseMatrix *C_stab = nullptr; // not owned
   #endif

   // data
   mfem::GridFunction *curr_vdot = nullptr; // not owned
   mfem::Array<int> *vel_ess_tdofs = nullptr; // not owned
   mfem::Array<int> *pres_ess_tdofs = nullptr; // not owned

   mfem::FiniteElementSpace *vfes = nullptr; // not owned
   mfem::FiniteElementSpace *pfes = nullptr; // not owned
   // total fluid integrator
   bool total_integ = false;
   mfem::BlockNonlinearForm *FluidNF = nullptr; // not owned

   // options
   bool equal_order = false;

   // cv bcs
   bool cv_bcs = false;
   CVBC_Op *cvbc_op = nullptr; // not owned
   mfem::Vector *A3d_diag = nullptr; // not owned
   mfem::Vector *Abc_diag = nullptr; // not owned
   mfem::Vector *Atotal_diag = nullptr; // not owned
   mutable SumOp *sum_grad_op = nullptr; // owned

   // immersed fsi
   bool IFEM = false;
   #ifdef MFEM_USE_MPI
   mutable mfem::HypreParMatrix *Agrad = nullptr; // owned if IFEM
   #else
   mutable mfem::SparseMatrix *Agrad = nullptr; // owned if IFEM
   #endif
   immersed_solid::ImmersedSolid *im_solid = nullptr; // not owned
   mfem::Vector *f_fsi = nullptr; // not owned

   // direct solver
   bool direct_solver = false;
   #ifdef MFEM_USE_MPI
   mutable mfem::Array2D<const mfem::HypreParMatrix*> gradMat; // blocks not owned
   mutable mfem::HypreParMatrix *totalGradMat = nullptr; // owned
   #else
   mutable mfem::BlockMatrix *gradMat = nullptr; // owned
   mutable mfem::SparseMatrix *totalGradMat = nullptr; // owned
   #endif

   // monolithic preconditioner
   bool monolithic_prec = false;

   // timers
   mutable mfem::StopWatch *sw_fluidtan = nullptr; // not owned
   mutable mfem::StopWatch *sw_fluidres = nullptr; // not owned
public:
   // constructor
   ResidualOpFluid(mfem::Array<int> block_trueOffsets_, mfem::Array<int> &vel_tdofs_, mfem::Array<int> &pres_tdofs_,
                     bool equal_order_) : 
   Operator(block_trueOffsets_[2]), block_trueOffsets(block_trueOffsets_), equal_order(equal_order_), 
   vel_ess_tdofs(&vel_tdofs_), pres_ess_tdofs(&pres_tdofs_)
   {
      temp.SetSize(block_trueOffsets[2]);
      temp = 0.0;
   };

   void AddCVBC_Op(CVBC_Op *cv_, mfem::Vector *A3d_diag_, mfem::Vector *Abc_diag_, mfem::Vector *Atotal_diag_)
   {
      cvbc_op = cv_;
      A3d_diag = A3d_diag_;
      Abc_diag = Abc_diag_;
      Atotal_diag = Atotal_diag_;
      cv_bcs = true;
   }

   void SetTotalInteg(mfem::BlockNonlinearForm *NF, mfem::FiniteElementSpace *vfes_, mfem::FiniteElementSpace *pfes_)
   {
      total_integ = true;
      FluidNF = NF;
      vfes = vfes_;
      pfes = pfes_;
   };

   void SetVDotGF(mfem::GridFunction &vdot_gf_)
   {
      curr_vdot = &vdot_gf_;
   };

   void SetIFEM(ImmersedSolid &im_solid_, mfem::Vector &f_fsi_)
   {
      IFEM = true;
      im_solid = &im_solid_;
      f_fsi = &f_fsi_;
   };

   void UsingDirectSolver()
   {
      direct_solver = true;
      #ifdef MFEM_USE_MPI
      gradMat.SetSize(2,2);
      gradMat = nullptr;
      #endif
   }

   void UsingMonolithicPrec()
   {
      monolithic_prec = true;
      #ifdef MFEM_USE_MPI
      gradMat.SetSize(2,2);
      gradMat = nullptr;
      #endif
   }

   void SetTimers(mfem::StopWatch *sw_fluidtan_, mfem::StopWatch *sw_fluidres_)
   {
      sw_fluidtan = sw_fluidtan_;
      sw_fluidres = sw_fluidres_;
   }

   // Compute the residual
   virtual void Mult(const mfem::Vector &k, mfem::Vector &y) const;
   void Mult_PC(const mfem::Vector &y, const mfem::Vector &ydot, const mfem::Vector &ft, const mfem::Vector &fb, double am, mfem::Vector &res, bool line_search = false) const;

   // Compute the gradient
   virtual mfem::Operator &GetGradient(const mfem::Vector &k) const;
   mfem::Operator &GetGradient_PC(const mfem::Vector &y, const mfem::Vector &ydot, double dt, double am, double af) const;
   void ComputeGradient(const mfem::Vector &y) const;

   // destructor
   virtual ~ResidualOpFluid();
};

} // ins_solver