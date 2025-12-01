#include "st_venant_kirchoff.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

void StVK::EvalPK2Stress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   // Compute S = lamdba*tr(E)I + 2*mu*E
   // stress = Identity_Mat;
   // stress *= lambda->Eval(T,ip);
   // double trE = E.Trace();
   // stress *= trE;
   // stress.Add(2.0 * mu->Eval(T,ip),E);

   stress.Set(2.0 * mu_val,E);
   double trE = E.Trace();
   stress(0,0) += lambda_val * trE;
   stress(1,1) += lambda_val * trE;
   if (vdim == 3) stress(2,2) += lambda_val * trE;
}

void StVK::EvalCauchyStress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   // compute sigma = 1/J * (lam * trE - mu)*b + mu/J b^2
   // Mult(b,b,b2);
   // b2 *= mu_val / J;
   // stress = b2;
   // double trE = E.Trace();
   // double factor = (lambda_val * trE - mu_val) / J;
   // stress.Add(factor,b);

   Mult(b,b,b2);
   stress.Set(mu_val/J,b2);
   double factor = (lambda_val * E.Trace() - mu_val) / J;
   stress.Add(factor,b);
}

double StVK::EvalElasticityTensor(int i, int j, int k, int l)
{
   if (formulation_type == form_type::TotalLagrangian)
   {
      // evaluate material elasticity tensor
      return lambda_val*kron_delta(i,j)*kron_delta(k,l) + mu_val*( kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l)*kron_delta(j,k) );
   }
   else if (formulation_type == form_type::UpdatedLagrangian)
   {
      // evaluate spatial elasticity tensor
      return (lambda_val / J)*b(i,j)*b(k,l) + (2.0*mu_val / J)*(1.0 / 2.0)*( b(i,k)*b(j,l) + b(i,l)*b(j,k) );
   }
   else
   {
      MFEM_ABORT("Error: Unsupported Formulation Type");
      return 0.0;
   }
}

double StVK::ElasticityTensorNorm()
{
   MFEM_ASSERT(formulation_type == form_type::UpdatedLagrangian, "Error this method is only supported for UL approach");

   // this is a guess, not sure if this will work!

   return lambda_val + 2.0*mu_val;
}

void StVK::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   lambda_val = lambda->Eval(T,ip);
   mu_val = mu->Eval(T,ip);

   mat_consts.SetSize(2);
   mat_consts(0) = mu_val; // c0
   mat_consts(1) = lambda_val; // c1
}

void StVK::ResetMaterialConstants()
{
   mu_val = mat_consts(0);
   lambda_val = mat_consts(1);
}

