#include "nearly_incompressible_transverse_isotropic.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

void NearlyIncompressibleFiberReinforced::ComputeS_iso(DenseMatrix &stress) 
{
   stress.SetSize(vdim);

   ///////////////////
   // Neo-Hookean Part

   double fac = -c0_val * pow(detC,-1.0/d) * 1.0/d * trC;
   stress.Set(fac,Cinv);
   stress(0,0) += c0_val * pow(detC,-1.0/d);
   stress(1,1) += c0_val * pow(detC,-1.0/d);
   if (vdim == 3) stress(2,2) += c0_val * pow(detC,-1.0/d);

   ///////////////////
   // Fung Part
   double i4_bar = pow(J,-2.0/d) * C.InnerProduct(m,m) - 1.0;
   double fac2 = c1_val * pow(J,-2.0/d)*(exp(i4_bar) - 1.0);

   // outerproduct of fiber directions m x m
   stress(0,0) += fac2 * m(0)*m(0);
   stress(0,1) += fac2 * m(0)*m(1);
   stress(1,0) += fac2 * m(1)*m(0);
   stress(1,1) += fac2 * m(1)*m(1);

   if (vdim == 3)
   {
      stress(0,2) += fac2 * m(0)*m(2);
      stress(1,2) += fac2 * m(1)*m(2);

      stress(2,0) += fac2 * m(2)*m(0);
      stress(2,1) += fac2 * m(2)*m(1);
      stress(2,2) += fac2 * m(2)*m(2);
   }
}

void NearlyIncompressibleFiberReinforced::ComputeSigma_iso(DenseMatrix &stress)
{
   MFEM_ABORT("TO DO");

   stress.SetSize(vdim);
   
   stress.Set(c0_val * pow(J,-(2.0 + d)/d),b);
   stress(0,0) -= c0_val * pow(J,-(2.0 + d)/d) * 1.0/d * b.Trace();
   stress(1,1) -= c0_val * pow(J,-(2.0 + d)/d) * 1.0/d * b.Trace();
   if (vdim == 3) stress(2,2) -= c0_val * pow(J,-(2.0 + d)/d) * 1.0/d * b.Trace();
}

double NearlyIncompressibleFiberReinforced::ComputeCref_iso(int i, int j, int k, int l)
{
   // deviatoric part
   // Neo-Hookean Part
   double C_dev = 2.0*c0_val*pow(detC,-1.0/d) * ( 1.0/d * trC * I_4th - 1.0/d * kron_delta(i,j)*Cinv(k,l) - 1.0/d * Cinv(i,j)*kron_delta(k,l) + 1.0/d/d * trC*Cinv(i,j)*Cinv(k,l));

   // Fung Part
   double i4_bar = pow(J,-2.0/d) * C.InnerProduct(m,m) - 1.0;
   double fac2 = 2.0 * c1_val * pow(J,-2.0/d)*exp(i4_bar);
   C_dev += fac2 * ( pow(J,-2.0/d)* m(i)*m(j)*m(k)*m(l) - 1.0/d * m(i)*m(j)*Cinv(k,l) );

   return C_dev;
}

double NearlyIncompressibleFiberReinforced::ComputeC_iso(int i, int j, int k, int l)
{
   MFEM_ABORT("TO DO");

   // deviatoric part
   double C_dev = 2.0*c0_val*pow(J,-(2.0 + d)/d) * ( 1.0/d * b.Trace()*i_val - 1.0/d * b(i,j)*kron_delta(k,l) - 1.0/d*kron_delta(i,j)*b(k,l) + 1.0/d/d*b.Trace()*kron_delta(i,j)*kron_delta(k,l));

   return C_dev;
}

void NearlyIncompressibleFiberReinforced::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   c0_val = c0->Eval(T,ip);
   c1_val = c1->Eval(T,ip);
   kappa_val = kappa->Eval(T,ip);

   mat_consts.SetSize(3);
   mat_consts(0) = c0_val;
   mat_consts(1) = c1_val;
   mat_consts(2) = kappa_val;

   if (fiber_coeff)
   {
      fiber_coeff->Eval(m,T,ip);
   }
}

void NearlyIncompressibleFiberReinforced::ResetMaterialConstants()
{
   c0_val = mat_consts(0);
   c1_val = mat_consts(1);
   kappa_val = mat_consts(2);
}

