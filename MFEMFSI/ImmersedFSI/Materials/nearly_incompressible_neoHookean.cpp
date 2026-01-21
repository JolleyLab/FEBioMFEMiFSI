#include "nearly_incompressible_neoHookean.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

void NearlyIncompressibleNeoHookean::ComputeS_iso(DenseMatrix &stress) 
{
   stress.SetSize(vdim);

   double fac = -mu_val * pow(detC,-1.0/d) * 1.0/d * trC;
   stress.Set(fac,Cinv);
   stress(0,0) += mu_val * pow(detC,-1.0/d);
   stress(1,1) += mu_val * pow(detC,-1.0/d);
   if (vdim == 3) stress(2,2) += mu_val * pow(detC,-1.0/d);
}

void NearlyIncompressibleNeoHookean::ComputeSigma_iso(DenseMatrix &stress)
{
   stress.SetSize(vdim);
      
   stress.Set(mu_val * pow(J,-(2.0 + d)/d),b);
   stress(0,0) -= mu_val * pow(J,-(2.0 + d)/d) * 1.0/d * b.Trace();
   stress(1,1) -= mu_val * pow(J,-(2.0 + d)/d) * 1.0/d * b.Trace();
   if (vdim == 3) stress(2,2) -= mu_val * pow(J,-(2.0 + d)/d) * 1.0/d * b.Trace();
}

double NearlyIncompressibleNeoHookean::ComputeCref_iso(int i, int j, int k, int l)
{
   // deviatoric part
   double C_dev = 2.0*mu_val*pow(detC,-1.0/d) * ( 1.0/d * trC * I_4th - 1.0/d * kron_delta(i,j)*Cinv(k,l) - 1.0/d * Cinv(i,j)*kron_delta(k,l) + 1.0/d/d * trC*Cinv(i,j)*Cinv(k,l));

   return C_dev;
}

double NearlyIncompressibleNeoHookean::ComputeC_iso(int i, int j, int k, int l)
{
   // deviatoric part
   double C_dev = 2.0*mu_val*pow(J,-(2.0 + d)/d) * ( 1.0/d * b.Trace()*i_val - 1.0/d * b(i,j)*kron_delta(k,l) - 1.0/d*kron_delta(i,j)*b(k,l) + 1.0/d/d*b.Trace()*kron_delta(i,j)*kron_delta(k,l));

   return C_dev;
}

void NearlyIncompressibleNeoHookean::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   kappa_val = kappa->Eval(T,ip);
   mu_val = mu->Eval(T,ip);

   mat_consts.SetSize(2);
   mat_consts(0) = mu_val; // c0
   mat_consts(1) = kappa_val; // c1
}

void NearlyIncompressibleNeoHookean::ResetMaterialConstants()
{
   mu_val = mat_consts(0);
   kappa_val = mat_consts(1);
}

