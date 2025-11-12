#include "nearly_incompressible_annular_disk.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

NearlyIncompressibleAnnularDisk::NearlyIncompressibleAnnularDisk(int form_type, int vdim, int vol_type, Coefficient &kappa_, Coefficient &mu_, double &width_) : NearlyIncompressibleMaterial(form_type, vdim,vol_type, kappa_) 
{
   mu = &mu_;
   w = width_;

   transip.SetSize(vdim);
   transip = 0.0;

   Finv.SetSize(vdim);
   Finv = 0.0;

   Pstress.SetSize(vdim);
   Pstress = 0.0;

   ehat_theta.SetSize(vdim);
   ehat_theta = 0.0;

   et_et.SetSize(vdim);
   et_et = 0.0;

   temp_vec.SetSize(vdim);      
};

void NearlyIncompressibleAnnularDisk::ComputeS_iso(DenseMatrix &stress) 
{
   stress.SetSize(vdim);

   F.Mult(ehat_theta,temp_vec);

   Pstress(0,0) = mu_val*temp_vec(0)*ehat_theta(0);
   Pstress(0,1) = mu_val*temp_vec(0)*ehat_theta(1);
   Pstress(1,0) = mu_val*temp_vec(1)*ehat_theta(0);
   Pstress(1,1) = mu_val*temp_vec(1)*ehat_theta(1);

   // compute S = Finv*P;
   Mult(Finv,Pstress,stress);  
}

void NearlyIncompressibleAnnularDisk::ComputeSigma_iso(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   // compute P
   F.Mult(ehat_theta,temp_vec);

   Pstress(0,0) = mu_val*temp_vec(0)*ehat_theta(0);
   Pstress(0,1) = mu_val*temp_vec(0)*ehat_theta(1);
   Pstress(1,0) = mu_val*temp_vec(1)*ehat_theta(0);
   Pstress(1,1) = mu_val*temp_vec(1)*ehat_theta(1);

   // sigma = 1/J P*F^T
   Mult(Pstress,Ft,temp_stress);
   stress.Set(1.0/J,temp_stress);
}

double NearlyIncompressibleAnnularDisk::ComputeCref_iso(int i, int j, int k, int l)
{  
   return mu_val * kron_delta(i,k) * ehat_theta(l) * ehat_theta(j);
}

double NearlyIncompressibleAnnularDisk::ComputeC_iso(int i, int j, int k, int l)
{
   // F.Mult(ehat_theta,temp_vec); // dont recompute, integrator is setup to compute stress above then C
   return (mu_val / J ) * b(i,k) * temp_vec(j) * temp_vec(l);
}

void NearlyIncompressibleAnnularDisk::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   kappa_val = kappa->Eval(T,ip);
   mu_val = mu->Eval(T,ip);

   mat_consts.SetSize(2);
   mat_consts(0) = mu_val; // c0
   mat_consts(1) = kappa_val; // c1

   // compute Finv
   Finv = F;
   Finv.Invert();

   // compute theta for current integration point
   T.Transform(ip,transip);   
   // shift to center of disk (0.5,0.5)
   transip(0) -= 0.5;
   transip(1) -= 0.5;

   // compute ehat theta
   ehat_theta(0) = -transip(1)/transip.Norml2();
   ehat_theta(1) = transip(0)/transip.Norml2();

   // compute theta
   theta = atan2(transip(1),transip(0));
   // shift theta between 0 and 2pi
   if (theta < 0.0)
   {
      double theta_temp = -M_PI - theta;
      theta = -theta_temp + M_PI;
   }
}

void NearlyIncompressibleAnnularDisk::ResetMaterialConstants()
{
   mu_val = mat_consts(0);
   kappa_val = mat_consts(1);
}

