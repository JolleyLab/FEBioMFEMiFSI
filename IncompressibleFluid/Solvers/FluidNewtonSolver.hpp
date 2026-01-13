#pragma once

#include "mfem.hpp"
#include "../../IncompressibleFluid/ifs_bc_containers.hpp"
#include "../SpatialOperators/ResidualOpFluid.hpp"
#include "BlockOrthoSolver.hpp"
#include "../../ImmersedFSI/ImmersedSolid.hpp"

using namespace immersed_solid;

namespace ifs_solver
{
// predictor corrector newton solver for fluid mechanics
// class to extend MFEM newton solver with additional options
class FluidNewtonSolver : public mfem::NewtonSolver
{
protected:
   int min_iter = 0; // min number of iterations to do (0 by default)
   mfem::Array<int> block_trueOffsets;

   int my_rank = 0;

   // integrator
   bool total_integ = false;

   // time stepping parameters
   // gen alpha
   double gamma = 0.0;
   double beta = 0.0;
   double alpha_m = 0.0;
   double alpha_f = 0.0;

   // Solution vectors
   mfem::Vector *yn = nullptr;
   mfem::Vector *ynp1 = nullptr;
   mfem::Vector *ydotn = nullptr;
   mfem::Vector *ydotnp1 = nullptr;

   // midpoint values for gen alpha
   mfem::Vector *yg = nullptr;
   mfem::Vector *ydotg = nullptr;

   // Special Residual Operator Class for Fluid Mechanics
   ResidualOpFluid *res_op = nullptr;

   // weak bcs
   std::vector<VelNeumannBC_T> *vel_nbcs = nullptr;
   mfem::LinearForm *L_L_form = nullptr;
   std::vector<AccelTerm_T> *accel_terms = nullptr;
   mfem::LinearForm *f_form = nullptr;
   mfem::Vector *fb = nullptr;
   mfem::Vector *ft = nullptr;

   // vars
   mutable mfem::BlockVector rblock;
   mutable int it, final_it_ls_fs;
   mutable double norm0, norm, norm_goal, final_atol_ls_fs, final_rtol_ls_fs;
   mutable double norm0_mom, norm0_mass, norm_mom, norm_mass;
   mutable double rel_norm_mom, rel_norm_mass;
   mutable double rel_norm, rel_norm_ls;
   bool print_newton = true;

   // direct solver
   bool direct_solver = false;

   // cvbcs
   std::vector<WindkesselBC_T> *windkessel_bcs = nullptr;
   bool windkessel_bc_flag = false;
   // no need for constant R resistance bcs!
   bool resistance_bc_flag = false;

   // immersed fsi
   immersed_solid::ImmersedSolid *im_solid = nullptr;
   bool IFEM = false;
   mutable mfem::BlockVector yblock, delblock;
   bool no_error_nonconvergence = false;

   // timing
   mfem::StopWatch *sw_tan, *sw_res, *sw_linsolve, *sw_initiators, *sw_correctors, *sw_weakbcs;

public:
   FluidNewtonSolver(ResidualOpFluid *res_, mfem::Array<int> &block_trueOffsets_) : mfem::NewtonSolver(), res_op(res_), block_trueOffsets(block_trueOffsets_) {};

   #ifdef MFEM_USE_MPI
   FluidNewtonSolver(ResidualOpFluid *res_, mfem::Array<int> &block_trueOffsets_, MPI_Comm comm_) : NewtonSolver(comm_), res_op(res_), block_trueOffsets(block_trueOffsets_) {};
   #endif

   // extended mult method
   void Mult(const mfem::Vector &b, mfem::Vector &x) const override; // needed since virtual in base class
   void NewtonSolve(const double &dt, const double &time) const;

   // pass in timing
   void SetTimers(mfem::StopWatch *sw_tan_, mfem::StopWatch *sw_res_, mfem::StopWatch *sw_linsolve_, mfem::StopWatch *sw_initiators_, mfem::StopWatch *sw_correctors_, mfem::StopWatch *sw_weakbcs_)
   {
      sw_tan = sw_tan_;
      sw_res = sw_res_;
      sw_linsolve = sw_linsolve_;
      sw_initiators = sw_initiators_;
      sw_correctors = sw_correctors_;
      sw_weakbcs = sw_weakbcs_;
   };

   // Set weak bcs
   void WeakBCs(std::vector<VelNeumannBC_T> *vel_nbcs_, mfem::LinearForm *L_L_form_, mfem::Vector &ft_, std::vector<AccelTerm_T> *accel_terms_, mfem::LinearForm *f_form_, mfem::Vector &fb_)
   {
      vel_nbcs = vel_nbcs_;
      L_L_form = L_L_form_;
      accel_terms = accel_terms_;
      f_form = f_form_;

      ft = &ft_;
      fb = &fb_;

      #ifdef MFEM_USE_MPI
      my_rank = dynamic_cast<mfem::ParLinearForm *>(L_L_form)->ParFESpace()->GetParMesh()->GetMyRank();
      #else
      my_rank = 0;
      #endif
   };

   // add cv bcs
   void AddWindkesselBC(std::vector<WindkesselBC_T> *windkessel_bcs_)
   {
      windkessel_bcs = windkessel_bcs_;
      windkessel_bc_flag = true;
   };

   // solution vectors
   void SetSolutionVectors(mfem::Vector &yn_, mfem::Vector &ydotn_, mfem::Vector &ynp1_, mfem::Vector &ydotnp1_, mfem::Vector &yg_, mfem::Vector &ydotg_)
   {
      yn = &yn_;
      ydotn = &ydotn_;
      ynp1 = &ynp1_;
      ydotnp1 = &ydotnp1_;
      yg = &yg_;
      ydotg = &ydotg_;
   };

   // time stepping method
   void UseGenAlpha(double alpha_f_, double alpha_m_, double gamma_)
   {
      alpha_f = alpha_f_;
      alpha_m = alpha_m_;
      gamma = gamma_;
   };

   void UsingDirectSolver()
   {
      direct_solver = true;
   }

   void SetIFEM(ImmersedSolid &im_solid_)
   {
      IFEM = true;
      im_solid = &im_solid_;
   };

   void SetNoErrorNonConvergence()
   {
      no_error_nonconvergence = true;
   };

   void SetTotalInteg()
   {
      total_integ = true;
   }

   void PrintNewton(bool val)
   {
      print_newton = val;
   }

   // additional functions for predictor corrector implementation
   void Initiator() const;
   void SetWeakBCs(const double &time, double dt) const;
   void FormLHS(mfem::Vector &y, mfem::Vector &ydot, double dt) const;
   void FormRHS(mfem::Vector &y, mfem::Vector &ydot, double dt, mfem::Vector &res) const;
   void SolveLinearSystem() const;
   void Correctors(const double &dt, const double &time) const;
   bool CheckConvergence() const;
   
   void GetResiduals(double &norm_, double &rel_norm_, double &norm_mom_, double &rel_norm_mom_, double &norm_mass_, double &rel_norm_mass_) const
   {
      norm_ = norm;
      rel_norm_ = rel_norm;

      norm_mom_ = norm_mom;
      rel_norm_mom_ = rel_norm_mom;

      norm_mass_ = norm_mass;
      rel_norm_mass_ = rel_norm_mass;
   };
   void PrintNorms() const;

   void PrintNewtonSolveInfo() const;
   void PrintLinSolveInfo() const;

   void SetMinIter(int min_it)
   {
      min_iter = min_it;
   };

   ~FluidNewtonSolver() override
   {
      
   };
};

}