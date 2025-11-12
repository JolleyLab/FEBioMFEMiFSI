#pragma once

#include "mfem.hpp"

using namespace mfem;
using namespace std;

class RobinBCIntegrator : public NonlinearFormIntegrator
{
private:
   Vector shape, u, udot;
   DenseMatrix EF;

   double ks = 0.0;
   double cs = 0.0;

   ParGridFunction *curr_v = nullptr;

   double *curr_dt = nullptr;
   double alpha_f;
   double gamma;
   double beta;
   bool div_factor = false;

public:
   /// Constructs a boundary integrator with a given VectorCoefficient QG
   RobinBCIntegrator(ParGridFunction &curr_v_, double ks_, double cs_, double *curr_dt_, double alpha_f_, double gamma_, double beta_, bool div_factor_ = true) :
   curr_v(&curr_v_), curr_dt(curr_dt_), alpha_f(alpha_f_), gamma(gamma_), beta(beta_), div_factor(div_factor_), ks(ks_), cs(cs_)
   { };

   // CG terms
   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &trans,
                                      const Vector &elfun,
                                      Vector &elvect);

   virtual void AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &trans,
                                     const Vector &elfun,
                                     DenseMatrix &elmat);

   virtual ~RobinBCIntegrator() {};
};