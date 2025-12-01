#include "nearly_incompressible_material.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

void NearlyIncompressibleMaterial::EvalPK2Stress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   // Compute S = S_iso + S_vol
   ComputeS_iso(stress);
   ComputeS_vol(stress);
}

void NearlyIncompressibleMaterial::EvalCauchyStress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   // Compute sigma = sigma_iso + sigma_vol
   ComputeSigma_iso(stress);
   ComputeSigma_vol(stress);
}

double NearlyIncompressibleMaterial::EvalElasticityTensor(int i, int j, int k, int l)
{
   double val = 0.0;
   if (formulation_type == form_type::TotalLagrangian)
   {
      I_4th = (1.0 / 2.0)*( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) );

      // ISO
      val = ComputeCref_iso(i,j,k,l);

      // VOL
      val += ComputeCref_vol(i,j,k,l);

      return val;
   }
   else if (formulation_type == form_type::UpdatedLagrangian)
   {
      i_val = 0.5*(kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l) * kron_delta(j,k));

      // ISO
      val = ComputeC_iso(i,j,k,l);

      // VOL
      val += ComputeC_vol(i,j,k,l);

      return val;
   }
   else
   {
      MFEM_ABORT("Error: Unsupported Formulation Type");
      return 0.0;
   }
}

double NearlyIncompressibleMaterial::ElasticityTensorNorm()
{
   MFEM_ABORT("Not implemented");
   return 0.0;
}

void NearlyIncompressibleMaterial::SetAdditionalStrainMeasures(const Vector &u, const DenseMatrix& gradu)
{
   // Compute C=F^T*F
   Mult(Ft,F,C);

   detC = C.Det();
   trC = C.Trace();

   // DenseMatrixInverse Cinv_op;
   Cinv_op.SetOperator(C);
   Cinv_op.GetInverseMatrix(Cinv);

   SetNearlyIncStrainMeasures(u,gradu);
}

void NearlyIncompressibleMaterial::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   MFEM_ABORT("Not implemented for base class");
}

void NearlyIncompressibleMaterial::ResetMaterialConstants()
{
   MFEM_ABORT("Not implemented for base class");
}

// volumetric energy functions
void NearlyIncompressibleMaterial::ComputeS_vol(DenseMatrix &stress)
{
   if (vol_type == volumetric_energy::QUAD)
   {
      stress.Add(J * (J - 1.0) * kappa_val, Cinv);
   }
   else if (vol_type == volumetric_energy::ST91)
   {
      stress.Add(0.5*kappa_val*(J*J - 1.0), Cinv);
   }
   else if (vol_type == volumetric_energy::M94)
   {
      MFEM_ABORT("TO DO");
   }
   else
   {
      MFEM_ABORT("Undefined volumetric energy type");
   }
}

void NearlyIncompressibleMaterial::ComputeSigma_vol(DenseMatrix &stress)
{
   if (vol_type == volumetric_energy::QUAD)
   {
      stress.Add(kappa_val * (J - 1.0),Identity_Mat);
   }
   else if (vol_type == volumetric_energy::ST91)
   {
      stress.Add(0.5*kappa_val*(J*J - 1.0)/J, Identity_Mat);
   }
   else if (vol_type == volumetric_energy::M94)
   {
      MFEM_ABORT("TO DO");
   }
   else
   {
      MFEM_ABORT("Undefined volumetric energy type");
   }
}

double NearlyIncompressibleMaterial::ComputeCref_vol(int i, int j, int k, int l)
{
   if (vol_type == volumetric_energy::QUAD)
   {
      // pressure part
      double C_p = kappa_val*(J-1.0)*J*(Cinv(i,j)*Cinv(k,l) - 2.0*I_4th);

      // volumetric part
      double C_vol = kappa_val * J * J * (Cinv(i,j)*Cinv(k,l));

      return C_p + C_vol;
   }
   else if (vol_type == volumetric_energy::ST91)
   {
      return kappa_val*J*J*Cinv(i,j)*Cinv(k,l) - kappa_val*(J*J - 1.0)*I_4th;
   }
   else if (vol_type == volumetric_energy::M94)
   {
      MFEM_ABORT("TO DO");
      return 0.0;
   }
   else
   {
      MFEM_ABORT("Undefined volumetric energy type");
      return 0.0;
   }
}

double NearlyIncompressibleMaterial::ComputeC_vol(int i, int j, int k, int l)
{
   if (vol_type == volumetric_energy::QUAD)
   {
      // pressure part
      double C_p = kappa_val*(J-1.0)*(kron_delta(i,j)*kron_delta(k,l) - 2.0*i_val);

      // volumetric part
      double C_vol = kappa_val * J *(kron_delta(i,j)*kron_delta(k,l));
      
      return C_p + C_vol;
   }
   else if (vol_type == volumetric_energy::ST91)
   {
      return kappa_val*J*kron_delta(i,j)*kron_delta(k,l) - kappa_val*(J*J - 1.0)*i_val / J;
   }
   else if (vol_type == volumetric_energy::M94)
   {
      MFEM_ABORT("TO DO");
      return 0.0;
   }
   else
   {
      MFEM_ABORT("Undefined volumetric energy type");
      return 0.0;
   }
}
