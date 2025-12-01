#include "compressible_neoHookean.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

void CompressibleNeoHookean::EvalPK2Stress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   // Compute S = mu*(I - Cinv) + lambda*ln(J)*Cinv
   // stress = Identity_Mat;
   // stress -= Cinv;
   // stress *= mu_val;

   stress.Set(-mu_val,Cinv);
   stress(0,0) += mu_val;
   stress(1,1) += mu_val;
   if (vdim == 3) stress(2,2) += mu_val;

   stress.Add(lambda_val * log(J),Cinv);
}

void CompressibleNeoHookean::EvalCauchyStress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   // compute sigma = mu/J (b - I) + (lambda/J)(lnJ) I
   // stress = b;
   // stress -= Identity_Mat;
   // stress *= (mu_val / J);
   // stress.Add(factor,Identity_Mat);

   stress.Set(mu_val / J,b);
   stress(0,0) -= mu_val / J;
   stress(1,1) -= mu_val / J;
   if (vdim == 3) stress(2,2) -= mu_val / J;   

   double factor = (lambda_val * log(J)) / J;
   stress(0,0) += factor;
   stress(1,1) += factor;
   if (vdim == 3) stress(2,2) += factor;
}

double CompressibleNeoHookean::EvalElasticityTensor(int i, int j, int k, int l)
{
   if (formulation_type == form_type::TotalLagrangian)
   {
      // evaluate material elasticity tensor
      return lambda_val*Cinv(i,j)*Cinv(k,l) + 2.0*(mu_val - lambda_val*log(J))*(1.0 / 2.0)*( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) );
   }
   else if (formulation_type == form_type::UpdatedLagrangian)
   {
      double lambda_prime = lambda_val / J;
      double mu_prime = (mu_val - lambda_val * log(J)) / J;
      
      // evaluate spatial elasticity tensor
      return lambda_prime*kron_delta(i,j)*kron_delta(k,l) + mu_prime*( kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l)*kron_delta(j,k) );
   }
   else
   {
      MFEM_ABORT("Error: Unsupported Formulation Type");
      return 0.0;
   }
}

double CompressibleNeoHookean::ElasticityTensorNorm()
{
   MFEM_ASSERT(formulation_type == form_type::UpdatedLagrangian, "Error this method is only supported for UL approach");

   double lambda_prime = lambda_val / J;
   double mu_prime = (mu_val - lambda_val * log(J)) / J;

   return lambda_prime + 2.0*mu_prime;
}

void CompressibleNeoHookean::SetAdditionalStrainMeasures(const Vector &u, const DenseMatrix& gradu)
{
   // Compute C=F^T*F
   DenseMatrix Ft(vdim);
   Ft = F;
   Ft.Transpose();
   Mult(Ft,F,C);

   DenseMatrixInverse Cinv_op;
   Cinv_op.SetOperator(C);
   Cinv_op.GetInverseMatrix(Cinv);
}

void CompressibleNeoHookean::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   lambda_val = lambda->Eval(T,ip);
   mu_val = mu->Eval(T,ip);

   mat_consts.SetSize(2);
   mat_consts(0) = mu_val; // c0
   mat_consts(1) = lambda_val; // c1
}

void CompressibleNeoHookean::ResetMaterialConstants()
{
   mu_val = mat_consts(0);
   lambda_val = mat_consts(1);
}

