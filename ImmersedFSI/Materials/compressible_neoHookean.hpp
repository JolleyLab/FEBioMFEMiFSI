#pragma once

// #include "mfem.hpp"
#include "material_type.hpp"

using namespace std;
using namespace mat;

class CompressibleNeoHookean : public MaterialType
{
protected:
   // elastic coefficients
	mfem::Coefficient *lambda; // c1
	mfem::Coefficient *mu; // c0
   double mu_val, lambda_val;

   // additional strain measure
   mfem::DenseMatrix C, Cinv;

public:
   CompressibleNeoHookean(int form_type, int vdim, mfem::Coefficient &lambda_, mfem::Coefficient &mu_) : MaterialType(form_type, vdim)
   {
      lambda = &lambda_;
      mu = &mu_;

      C.SetSize(vdim);
      C = 0.0;

      Cinv.SetSize(vdim);
      Cinv = 0.0;
   };

   // pk2 at quad pt
   virtual void EvalPK2Stress(mfem::DenseMatrix &stress) override;

   // cauchy stress at quad pt
   virtual void EvalCauchyStress(mfem::DenseMatrix &stress) override;

   // method to evaluate Elasticity Tensor
   virtual double EvalElasticityTensor(int i, int j, int k, int l) override;

   // virtual function called within SetStrainMeasures to compute additional quantities needed specific to a material model
   virtual void SetAdditionalStrainMeasures(const mfem::Vector &u, const mfem::DenseMatrix& gradu) override;

   // method to compute l2 norm of elasticity tensor
   virtual double ElasticityTensorNorm() override;

   virtual void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T) override;
   virtual void ResetMaterialConstants() override;
   virtual double GetKappa() override
   {
      return lambda_val + 2.0*mu_val/3.0;   
   };

   ~CompressibleNeoHookean() {};
};