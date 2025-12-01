#include "nearly_incompressible_Lee-Sacks.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

void NearlyIncompressibleLeeSacks::ComputeS_iso(DenseMatrix &stress) 
{
   stress.SetSize(vdim);

   // (A) (Neo-Hookean Part)
   double fac = -c0_val * pow(detC,-1.0/d) * 1.0/d * trC;
   stress.Set(fac,Cinv);
   stress(0,0) += c0_val * pow(detC,-1.0/d);
   stress(1,1) += c0_val * pow(detC,-1.0/d);
   if (vdim == 3) stress(2,2) += c0_val * pow(detC,-1.0/d);

   // (B) (Isotropic Part)
   tr_val = pow(J,-2.0/d)*trC - 3.0;
   double fac2 = 2.0*c1_val*delta_val*c2_val*tr_val*exp(c2_val*tr_val*tr_val) * pow(J,-2.0/d);
   stress(0,0) += fac2;
   stress(1,1) += fac2;
   if (vdim == 3) stress(2,2) += fac2;

   // (C) (Anisotropic Part)
   i4_bar = pow(J,-2.0/d) * C.InnerProduct(m,m) - 1.0;
   double fac3 = c1_val*(1.0 - delta_val)*exp(c3_val * i4_bar * i4_bar)*2.0*c3_val*i4_bar*pow(J,-2.0/d);

   // outerproduct of fiber directions m x m
   stress(0,0) += fac3 * m(0)*m(0);
   stress(0,1) += fac3 * m(0)*m(1);
   stress(1,0) += fac3 * m(1)*m(0);
   stress(1,1) += fac3 * m(1)*m(1);

   if (vdim == 3)
   {
      stress(0,2) += fac3 * m(0)*m(2);
      stress(1,2) += fac3 * m(1)*m(2);

      stress(2,0) += fac3 * m(2)*m(0);
      stress(2,1) += fac3 * m(2)*m(1);
      stress(2,2) += fac3 * m(2)*m(2);
   }
}

void NearlyIncompressibleLeeSacks::ComputeSigma_iso(DenseMatrix &stress)
{
   MFEM_ABORT("TO DO");
}

double NearlyIncompressibleLeeSacks::ComputeCref_iso(int i, int j, int k, int l)
{
   // (A) Neo-Hookean Part
   double C_dev = 2.0*c0_val*pow(detC,-1.0/d) * ( 1.0/d * trC * I_4th - 1.0/d * kron_delta(i,j)*Cinv(k,l) - 1.0/d * Cinv(i,j)*kron_delta(k,l) + 1.0/d/d * trC*Cinv(i,j)*Cinv(k,l));

   // (B) Isotropic Part
   double fac0 = 4.0*c1_val*delta_val*c2_val*kron_delta(i,j);
   double t1 = pow(J,-4.0/d)*exp(c2_val*tr_val*tr_val) * (1.0 + 2.0*c2_val*tr_val*tr_val)*kron_delta(k,l);
   double t2 = -1.0/d*pow(J,-2.0/d)*tr_val*exp(c2_val*tr_val*tr_val)*Cinv(k,l);

   C_dev += fac0 * (t1 + t2); 

   // (C) Anisotropic Part
   double fac1 = 4.0*c1_val*(1.0-delta_val)*c3_val*m(i)*m(j);
   double a1 = pow(J,-4.0/d)*exp(c3_val * i4_bar * i4_bar) * (1.0 + 2.0*c3_val*i4_bar*i4_bar) * m(k) * m(l);
   double a2 = -1.0/d*pow(J,-2.0/d)*i4_bar*exp(c3_val*i4_bar*i4_bar)*Cinv(k,l);

   C_dev += fac1 * (a1 + a2);

   return C_dev;
}

double NearlyIncompressibleLeeSacks::ComputeC_iso(int i, int j, int k, int l)
{
   // // deviatoric part
   // double C_dev = 2.0*mu_val*pow(J,-(2.0 + d)/d) * ( 1.0/d * b.Trace()*i_val - 1.0/d * b(i,j)*kron_delta(k,l) - 1.0/d*kron_delta(i,j)*b(k,l) + 1.0/d/d*b.Trace()*kron_delta(i,j)*kron_delta(k,l));

   // return C_dev;

   MFEM_ABORT("TO DO");
   return 0.0;
}

void NearlyIncompressibleLeeSacks::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   c0_val = c0->Eval(T,ip);
   c1_val = c1->Eval(T,ip);
   c2_val = c2->Eval(T,ip);
   c3_val = c3->Eval(T,ip);
   kappa_val = kappa->Eval(T,ip);

   mat_consts.SetSize(5);
   mat_consts(0) = c0_val;
   mat_consts(1) = c1_val;
   mat_consts(2) = c2_val;
   mat_consts(3) = c3_val;
   mat_consts(4) = kappa_val;

   if (fiber_coeff)
   {
      fiber_coeff->Eval(m,T,ip);
   }
}

void NearlyIncompressibleLeeSacks::ResetMaterialConstants()
{
   c0_val = mat_consts(0);
   c1_val = mat_consts(1);
   c2_val = mat_consts(2);
   c3_val = mat_consts(3);
   kappa_val = mat_consts(4);
}

