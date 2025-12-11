#pragma once

#include "mfem.hpp"

using namespace std;

// Nonlinear integrator for backflow stabilization term, -1.0*beta*rho*<v,(u*n)_ u>
class BackflowStabIntegrator : public mfem::NonlinearFormIntegrator
{
private:
   double beta = 0.2; // default
   mfem::Coefficient *coeff;
   mfem::VectorCoefficient *vec_coeff;
   bool ale;

   mfem::Vector shape, nor, u, uhat, un_u, u_bf;
   mfem::DenseMatrix EF, a11;
public:

   BackflowStabIntegrator() : coeff(NULL), vec_coeff(nullptr), ale(false) {};
   BackflowStabIntegrator(mfem::Coefficient &coeff_) : coeff(&coeff_), vec_coeff(nullptr), ale(false) {};
   BackflowStabIntegrator(mfem::Coefficient &coeff_, double beta_) : coeff(&coeff_), vec_coeff(nullptr), ale(false), beta(beta_) {};
   BackflowStabIntegrator(mfem::VectorCoefficient &vcoeff_, bool ale_) : coeff(nullptr), vec_coeff(&vcoeff_), ale(ale_) {};
   BackflowStabIntegrator(mfem::VectorCoefficient &vec_coeff_, mfem::Coefficient &coeff_, bool ale_) : coeff(&coeff_), vec_coeff(&vec_coeff_), ale(ale_) {};
   BackflowStabIntegrator(mfem::VectorCoefficient &vec_coeff_, mfem::Coefficient &coeff_, bool ale_, double beta_) : coeff(&coeff_), vec_coeff(&vec_coeff_), ale(ale_), beta(beta_) {};


   // CG terms
   virtual void AssembleElementVector(const mfem::FiniteElement &el,
                                      mfem::ElementTransformation &trans,
                                      const mfem::Vector &elfun,
                                      mfem::Vector &elvect);

   virtual void AssembleElementGrad(const mfem::FiniteElement &el,
                                     mfem::ElementTransformation &trans,
                                     const mfem::Vector &elfun,
                                     mfem::DenseMatrix &elmat);
};

class BlockBackflowStabIntegrator : public mfem::BlockNonlinearFormIntegrator
{
private:
   double beta = 0.2; // default
   mfem::Coefficient *coeff = nullptr;
   mfem::VectorGridFunctionCoefficient *mesh_velo = nullptr;
   bool ale;

   mfem::Vector shape, nor, u, uhat, un_u, u_bf;
   mfem::DenseMatrix EF;

   double *dt = nullptr;
   double alpha_f, gamma;
   bool gen_alpha = false;
public:

   BlockBackflowStabIntegrator(double *dt_, double alpha_f_, double gamma_, mfem::VectorGridFunctionCoefficient &vec_coeff_, mfem::Coefficient &coeff_, bool ale_, double beta_, bool gen_alpha_)
   : coeff(&coeff_), mesh_velo(&vec_coeff_), ale(ale_), beta(beta_), dt(dt_), alpha_f(alpha_f_), gamma(gamma_), gen_alpha(gen_alpha_) {};

   BlockBackflowStabIntegrator(double *dt_, double alpha_f_, double gamma_, mfem::Coefficient &coeff_, double beta_, bool gen_alpha_)
   : coeff(&coeff_), mesh_velo(nullptr), ale(false), beta(beta_), dt(dt_), alpha_f(alpha_f_), gamma(gamma_), gen_alpha(gen_alpha_) {};


   // CG terms
   virtual void AssembleElementVector(const mfem::Array<const mfem::FiniteElement *> &el,
                                       mfem::ElementTransformation &Tr,
                                       const mfem::Array<const mfem::Vector *> &elfun,
                                       const mfem::Array<mfem::Vector *> &elvec);

   virtual void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*> &el,
                                     mfem::ElementTransformation &Tr,
                                     const mfem::Array<const mfem::Vector *> &elfun,
                                     const mfem::Array2D<mfem::DenseMatrix *> &elmats);
};