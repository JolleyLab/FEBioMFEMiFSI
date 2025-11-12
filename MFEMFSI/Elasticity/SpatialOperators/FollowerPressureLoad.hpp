#pragma once

#include "mfem.hpp"

using namespace mfem;
using namespace std;

// Nonlinear integrator for Follower Pressure Load, assumes applied pressure is independent of the displacement 
// Currently implemented for 3D only!! 2D TO DO
class FollowerPressureLoad : public NonlinearFormIntegrator
{
private:
   Coefficient *pcoeff; // scalar pressure
   Coefficient *coeff = nullptr;

   double *curr_dt = nullptr;
   double alpha_f, beta;

   // for TL
   // for building current parent->phys element transformation
   IsoparametricTransformation *curr_trans = nullptr;
   int formulation_type = form_type::UpdatedLagrangian;
   DenseMatrix EF;

public:
   FollowerPressureLoad(Coefficient *coeff_) : pcoeff(coeff_), coeff(nullptr) {};
   FollowerPressureLoad(Coefficient *coeff_, Coefficient *coeff_const_) : pcoeff(coeff_), coeff(coeff_const_) {};
   FollowerPressureLoad(int form_, Coefficient *coeff_) : pcoeff(coeff_), coeff(nullptr), formulation_type(form_) {};
   FollowerPressureLoad(int form_, Coefficient *coeff_, Coefficient *coeff_const_) : pcoeff(coeff_), coeff(coeff_const_), formulation_type(form_) {};
   FollowerPressureLoad(int form_, Coefficient *coeff_, Coefficient *coeff_const_, double *curr_dt_, double alpha_f_, double beta_) 
   : pcoeff(coeff_), coeff(coeff_const_), formulation_type(form_), curr_dt(curr_dt_), alpha_f(alpha_f_), beta(beta_) {};

   // CG
   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &trans,
                                      const Vector &elfun,
                                      Vector &elvect);

   virtual void AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &trans,
                                     const Vector &elfun,
                                     DenseMatrix &elmat);

   // DG
   virtual void AssembleFaceVector(const FiniteElement &el1, const FiniteElement &el2,
                                   FaceElementTransformations &Tr, const Vector &elfun, Vector &elvect);

   virtual void AssembleFaceGrad(const FiniteElement &el1, const FiniteElement &el2,
                                   FaceElementTransformations &Tr, const Vector &elfun, DenseMatrix &elmat);

   // Building current transformation
   void BuildCurrentTransformation(ElementTransformation &T, const IntegrationPoint &ip, DenseMatrix *disp_pm);

   enum form_type {TotalLagrangian = 0, UpdatedLagrangian = 1};
};