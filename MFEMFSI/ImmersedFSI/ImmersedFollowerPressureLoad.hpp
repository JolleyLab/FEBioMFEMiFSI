#pragma once

#include "mfem.hpp"

using namespace std;

// Nonlinear integrator for Follower Pressure Load, assumes applied pressure is independent of the displacement 
// Currently implemented for 3D only!! 2D TO DO
class ImmersedFollowerPressureLoad : public mfem::NonlinearFormIntegrator
{
private:
   mfem::Coefficient *pcoeff; // scalar pressure
   mfem::Coefficient *coeff = nullptr;

   double *dt = nullptr;
   double alpha_f, beta;

   // for TL
   // for building current parent->phys element transformation
   mfem::IsoparametricTransformation *curr_trans = nullptr;
   int formulation_type = form_type::UpdatedLagrangian;
   mfem::DenseMatrix EF;

public:
   ImmersedFollowerPressureLoad(int form_, mfem::Coefficient *coeff_, mfem::Coefficient *coeff_const_, double *curr_dt_, double alpha_f_, double beta_) 
   : pcoeff(coeff_), coeff(coeff_const_), formulation_type(form_), dt(curr_dt_), alpha_f(alpha_f_), beta(beta_) {};

   virtual void AssembleElementVector(const mfem::FiniteElement &el,
                                      mfem::ElementTransformation &trans,
                                      const mfem::Vector &elfun,
                                      mfem::Vector &elvect);

   virtual void AssembleElementGrad(const mfem::FiniteElement &el,
                                     mfem::ElementTransformation &trans,
                                     const mfem::Vector &elfun,
                                     mfem::DenseMatrix &elmat);

   // Building current transformation
   void BuildCurrentTransformation(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip, mfem::DenseMatrix *disp_pm);

   enum form_type {TotalLagrangian = 0, UpdatedLagrangian = 1};
};