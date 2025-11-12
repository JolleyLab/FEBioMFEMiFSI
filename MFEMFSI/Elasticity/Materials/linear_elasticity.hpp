#pragma once

#include "material_type.hpp"

using namespace std;
using namespace mat;

class LinearElasticity : public MaterialType
{
protected:
   // elastic coefficients
	mfem::Coefficient *lambda; // c1
	mfem::Coefficient *mu; // c0

   double lambda_val, mu_val;

   // infinitesimal strain
   mfem::DenseMatrix e;
public:
   LinearElasticity(int form_type, int vdim, mfem::Coefficient &lambda_, mfem::Coefficient &mu_) : MaterialType(form_type, vdim)
   {
      lambda = &lambda_;
      mu = &mu_;

      e.SetSize(vdim);
      e = 0.0;
   };

   // method to evaluate PK2 Stress = Cauchy Stress for LinearElasticity
   virtual void EvalPK2Stress(mfem::DenseMatrix &stress) override;

   // method to evaluate PK1 Stress = Cauchy Stress for LinearElasticity
   virtual void EvalPK1Stress(mfem::DenseMatrix &stress) override;

   // method to evaluate Cauchy Stress given PK1 or directly evaluate Cauchy stress
   virtual void EvalCauchyStress(mfem::DenseMatrix &stress) override;

   // method to evaluate Elasticity Tensor
   virtual double EvalElasticityTensor(int i, int j, int k, int l) override;

   virtual void SetAdditionalStrainMeasures(const mfem::Vector &u, const mfem::DenseMatrix& gradu) override
   {
      // e = sym(gradu)
      e = gradu;
	  mfem::DenseMatrix graduT(vdim);
      graduT = gradu;
      graduT.Transpose();  
      e += graduT;    
      e *= 0.5;
   }

   // method to compute l2 norm of elasticity tensor
   virtual double ElasticityTensorNorm() override;

   virtual void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T) override;
   virtual void ResetMaterialConstants() override;

   virtual double GetKappa() override
   {
      return lambda_val + 2.0*mu_val/3.0;   
   };
   
   // // assembles local element vector for stress term
   // virtual void AssembleElementVector(const FiniteElement &el,
   //                                     ElementTransformation &trans,
   //                                     const Vector &elfun,
   //                                     Vector &elvect) override;

   // void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat);

   // // assembles local element grad matrix for stress term
   // virtual void AssembleElementGrad(const FiniteElement &el,
   //                                     ElementTransformation &trans,
   //                                     const Vector &elfun,
   //                                     DenseMatrix &elmat) override;
};