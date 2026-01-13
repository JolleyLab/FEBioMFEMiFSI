 #pragma once

// #include "mfem.hpp"
#include "material_type.hpp"
#include "nearly_incompressible_material.hpp"

using namespace std;
using namespace mat;

// nearly incompressible neo hookean model with fung type fiber reinforcement (transversely isotropic)
// see https://doi.org/10.1016/0045-7825(96)01035-3 
class NearlyIncompressibleFiberReinforced : public NearlyIncompressibleMaterial
{
protected:
   // elastic coefficients
	mfem::Coefficient *c1 = nullptr; // c1
	mfem::Coefficient *c0 = nullptr; // c0

	mfem::Vector m; // fiber direction
	mfem::VectorCoefficient *fiber_coeff = nullptr;

   double c0_val = 0.0;
   double c1_val = 0.0;

public:
   NearlyIncompressibleFiberReinforced(int form_type, int vdim, int vol_type, mfem::Coefficient &kappa_, mfem::Coefficient &c0_, mfem::Coefficient &c1_, mfem::Vector &fiber_dir) : NearlyIncompressibleMaterial(form_type, vdim,vol_type, kappa_)
   {
      c1 = &c1_;
      c0 = &c0_;

      // check if unit vector if not normalize it
      double norm_fiber = fiber_dir.Norml2();
      if (abs(norm_fiber - 1.0) < 1e-12) m = fiber_dir;
      else
      {
         m = fiber_dir;
         m /= norm_fiber;
      }

      m = fiber_dir;
   };
   NearlyIncompressibleFiberReinforced(int form_type, int vdim, int vol_type, mfem::Coefficient &kappa_, mfem::Coefficient &c0_, mfem::Coefficient &c1_, mfem::VectorCoefficient &fiber_dir) : NearlyIncompressibleMaterial(form_type, vdim,vol_type, kappa_)
   {
      c1 = &c1_;
      c0 = &c0_;

      // assume unit fiber direction
      fiber_coeff = &fiber_dir;
   };

   virtual void ComputeS_iso(mfem::DenseMatrix &stress) override;
   virtual void ComputeSigma_iso(mfem::DenseMatrix &stress) override;

   virtual double ComputeCref_iso(int i, int j, int k, int l) override;
   virtual double ComputeC_iso(int i, int j, int k, int l) override;

   virtual void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T) override;
   virtual void ResetMaterialConstants() override;

   ~NearlyIncompressibleFiberReinforced() {};
};