 #pragma once

// #include "mfem.hpp"
#include "material_type.hpp"
#include "nearly_incompressible_material.hpp"

using namespace std;
using namespace mat;

// nearly incompressible neo hookean model where incompressibility is enforced via an augmented strain energy
class NearlyIncompressibleNeoHookean : public NearlyIncompressibleMaterial
{
protected:
   // elastic coefficients
	mfem::Coefficient *mu; // c0

   double mu_val;

public:
   NearlyIncompressibleNeoHookean(int form_type, int vdim, int vol_type, mfem::Coefficient &kappa_, mfem::Coefficient &mu_) : NearlyIncompressibleMaterial(form_type, vdim,vol_type, kappa_)
   {
      mu = &mu_;
   };

   virtual void ComputeS_iso(mfem::DenseMatrix &stress) override;
   virtual void ComputeSigma_iso(mfem::DenseMatrix &stress) override;

   virtual double ComputeCref_iso(int i, int j, int k, int l) override;
   virtual double ComputeC_iso(int i, int j, int k, int l) override;

   virtual void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T) override;
   virtual void ResetMaterialConstants() override;

   ~NearlyIncompressibleNeoHookean() {};
};