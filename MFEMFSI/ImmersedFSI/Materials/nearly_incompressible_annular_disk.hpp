#pragma once

// #include "mfem.hpp"
#include "nearly_incompressible_material.hpp"

using namespace std;
using namespace mat;

// nearly incompressible annular disk for simple test problem
class NearlyIncompressibleAnnularDisk : public NearlyIncompressibleMaterial
{
protected:
   // elastic coefficients
	mfem::Coefficient *mu; // c0
   double mu_val, w;
   mfem::Vector transip, ehat_theta, temp_vec;
   double theta;
   mfem::DenseMatrix Finv, Pstress, et_et, temp;
public:
   NearlyIncompressibleAnnularDisk(int form_type, int vdim, int vol_type, mfem::Coefficient &kappa_, mfem::Coefficient &mu_, double &width_);

   virtual void ComputeS_iso(mfem::DenseMatrix &stress) override;
   virtual void ComputeSigma_iso(mfem::DenseMatrix &stress) override;

   virtual double ComputeCref_iso(int i, int j, int k, int l) override;
   virtual double ComputeC_iso(int i, int j, int k, int l) override;

   virtual void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T) override;
   virtual void ResetMaterialConstants() override;

   double Q_matrix(int i, int j)
   {
      if (i == 0 && j == 0)
      {
         return cos(theta);
      }
      else if (i == 0 && j == 1)
      {
         return -sin(theta);
      }
      else if (i == 1 && j == 0)
      {
         return sin(theta);
      }
      else if (i == 1 && j == 1)
      {
         return cos(theta);
      }
      else
      {
         MFEM_ABORT("Error: invalid matrix entry");
      }
   };

   ~NearlyIncompressibleAnnularDisk() {};
};