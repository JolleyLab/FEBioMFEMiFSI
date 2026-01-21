#include "linear_elasticity.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

void LinearElasticity::EvalPK2Stress(DenseMatrix &stress)
{
   EvalCauchyStress(stress);
}

void LinearElasticity::EvalPK1Stress(DenseMatrix &stress)
{
   EvalCauchyStress(stress);
}

void LinearElasticity::EvalCauchyStress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   stress.Set(e.Trace()*lambda_val,Identity_Mat);
   stress.Add(2.0 * mu_val,e);
}

double LinearElasticity::EvalElasticityTensor(int i, int j, int k, int l)
{
   if (formulation_type == form_type::TotalLagrangian)
   {
      // evaluate material elasticity tensor
      return lambda_val*kron_delta(i,j)*kron_delta(k,l) + mu_val*( kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l)*kron_delta(j,k) );

   }
   else if (formulation_type == form_type::UpdatedLagrangian)
   {
      // evaluate spatial elasticity tensor
      return lambda_val*kron_delta(i,j)*kron_delta(k,l) + mu_val*( kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l)*kron_delta(j,k) );
   }
   else
   {
      MFEM_ABORT("Error: Unsupported Formulation Type");
      return 0.0;
   }
}

double LinearElasticity::ElasticityTensorNorm()
{
   MFEM_ASSERT(formulation_type == form_type::UpdatedLagrangian, "Error this method is only supported for UL approach");

   return lambda_val + 2.0*mu_val;
}

void LinearElasticity::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   lambda_val = lambda->Eval(T,ip);
   mu_val = mu->Eval(T,ip);

   mat_consts.SetSize(2);
   mat_consts(0) = mu_val; // c0
   mat_consts(1) = lambda_val; // c1
}

void LinearElasticity::ResetMaterialConstants()
{
   mu_val = mat_consts(0);
   lambda_val = mat_consts(1);
}
                                 