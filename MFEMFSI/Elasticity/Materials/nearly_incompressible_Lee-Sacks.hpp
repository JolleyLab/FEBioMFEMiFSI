 #pragma once

// #include "mfem.hpp"
#include "material_type.hpp"
#include "nearly_incompressible_material.hpp"

using namespace std;
using namespace mat;

// nearly incompressible neo hookean model where incompressibility is enforced via an augmented strain energy
class NearlyIncompressibleLeeSacks : public NearlyIncompressibleMaterial
{
protected:
   // elastic coefficients
   mfem::Coefficient *c0;
   mfem::Coefficient *c1;
   mfem::Coefficient *c2;
   mfem::Coefficient *c3;

   double c0_val, c1_val, c2_val, c3_val;

   // fiber direction (in global specimen coordinates, NOT material coordinates)
   // UNIT VECTOR
   mfem::Vector m;
   mfem::VectorCoefficient *fiber_coeff = nullptr;

   // parameter
   double delta_val;

   // temps
   double tr_val;
   double i4_bar;

public:
   NearlyIncompressibleLeeSacks(int form_type, int vdim, int vol_type, mfem::Coefficient &kappa_, mfem::Coefficient &C0_, mfem::Coefficient &C1_, mfem::Coefficient &C2_, mfem::Coefficient &C3_, mfem::Vector &fiber_dir, double &delta_val_) : NearlyIncompressibleMaterial(form_type, vdim,vol_type, kappa_)
   {
      c0 = &C0_;
      c1 = &C1_;
      c2 = &C2_;
      c3 = &C3_;

      // check if unit vector if not normalize it
      double norm_fiber = fiber_dir.Norml2();
      if (abs(norm_fiber - 1.0) < 1e-12) m = fiber_dir;
      else
      {
         m = fiber_dir;
         m /= norm_fiber;
      }

      m = fiber_dir;


      delta_val = delta_val_;
   };
   NearlyIncompressibleLeeSacks(int form_type, int vdim, int vol_type, mfem::Coefficient &kappa_, mfem::Coefficient &C0_, mfem::Coefficient &C1_, mfem::Coefficient &C2_, mfem::Coefficient &C3_, mfem::VectorCoefficient &fiber_dir, double &delta_val_) : NearlyIncompressibleMaterial(form_type, vdim,vol_type, kappa_)
   {
      c0 = &C0_;
      c1 = &C1_;
      c2 = &C2_;
      c3 = &C3_;

      // check if unit vector if not normalize it
      fiber_coeff = &fiber_dir;

      delta_val = delta_val_;
   };

   virtual void ComputeS_iso(mfem::DenseMatrix &stress) override;
   virtual void ComputeSigma_iso(mfem::DenseMatrix &stress) override;

   virtual double ComputeCref_iso(int i, int j, int k, int l) override;
   virtual double ComputeC_iso(int i, int j, int k, int l) override;

   virtual void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T) override;
   virtual void ResetMaterialConstants() override;

   ~NearlyIncompressibleLeeSacks() {};
};