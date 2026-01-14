 #pragma once

// #include "mfem.hpp"
#include "material_type.hpp"

using namespace std;
using namespace mat;

// nearly incompressible material model
/*
   Based on additive decomp of strain energy into volumetric and isochoric parts
   Base class implements multiple volumetric energy functions
   Derived classes will define isochoric parts
*/
class NearlyIncompressibleMaterial : public MaterialType
{
protected:
   // volumetric energy type
   int vol_type;

   // additional strain measure
   mfem::DenseMatrix C, Cinv;
   mfem::DenseMatrixInverse Cinv_op;
   double detC, trC;

   mfem::Coefficient *kappa = nullptr;

   double kappa_val = 0.0;
   double I_4th = 0.0;
   double i_val = 0.0;

public:
   NearlyIncompressibleMaterial(int form_type, int vdim) : MaterialType(form_type, vdim) 
   {
      C.SetSize(vdim);
      C = 0.0;

      Cinv.SetSize(vdim);
      Cinv = 0.0;
   };
   NearlyIncompressibleMaterial(int form_type, int vdim, int vol_type_, mfem::Coefficient &kappa_) : MaterialType(form_type, vdim), vol_type(vol_type_)
   {
      C.SetSize(vdim);
      C = 0.0;

      Cinv.SetSize(vdim);
      Cinv = 0.0;

      kappa = &kappa_;
   };   

   void SetKappaCoefficient(mfem::Coefficient &kappa_)
   {
      kappa = &kappa_;
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

   virtual double GetKappa() override
   {
      return kappa_val;   
   };

   ////////////////////////////////////////////////////
   // functions that need to be overloaded
   virtual void SetNearlyIncStrainMeasures(const mfem::Vector &u, const mfem::DenseMatrix& gradu)
   {
      // implement as needed for derived classes
      // called within SetAdditionalStrainMeasures
   };

   virtual void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T) override;
   virtual void ResetMaterialConstants() override;

   // functions to compute isochoric part
   // overloaded in derived classes
   virtual void ComputeS_iso(mfem::DenseMatrix &stress)
   {
      MFEM_ABORT("Not implemented for base class");
   };

   virtual void ComputeSigma_iso(mfem::DenseMatrix &stress)
   {
      MFEM_ABORT("Not implemented for base class");
   };

   virtual double ComputeCref_iso(int i, int j, int k, int l)
   {
      MFEM_ABORT("Not implemented for base class");
   };

   virtual double ComputeC_iso(int i, int j, int k, int l)
   {
      MFEM_ABORT("Not implemented for base class");
   };

   ////////////////////////////////////////////////////
   // volumetric energy options and functions
   enum volumetric_energy {QUAD = 0, ST91 = 1, M94 = 2};

   void ComputeS_vol(mfem::DenseMatrix &stress);
   void ComputeSigma_vol(mfem::DenseMatrix &stress);
   double ComputeCref_vol(int i, int j, int k, int l);
   double ComputeC_vol(int i, int j, int k, int l);


   virtual ~NearlyIncompressibleMaterial() {};
};