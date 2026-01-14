#include "material_type.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

// method to evaluate PK2 Stress
void MaterialType::EvalPK2Stress(DenseMatrix &stress)
{
   MFEM_ABORT("MaterialType: Not implemented for base class");
}

// method to evaluate PK1 Stress
void MaterialType::EvalPK1Stress(DenseMatrix &stress)
{
   // compute S
   EvalPK2Stress(stress);

   // P = FS
   Mult(F,stress,temp_stress);
   stress = temp_stress;   
}

// method to evaluate Cauchy Stress
void MaterialType::EvalCauchyStress(DenseMatrix &stress)
{
   // compute P
   EvalPK1Stress(stress);

   // sigma = 1/J P*F^T
   Mult(stress,Ft,temp_stress);
   stress.Set(1.0/J,temp_stress);
}

double MaterialType::kron_delta(int i, int j)
{
   if (i == j)
   {
      return 1.0;
   }
   else
   {
      return 0.0;
   }
}

double MaterialType::EvalElasticityTensor(int i, int j, int k, int l)
{
   MFEM_ABORT("MaterialType: Not implemented for base class");
   return 0.0;
}

double MaterialType::ElasticityTensorNorm()
{
   MFEM_ABORT("MaterialType: Not implemented for base class");
   return 0.0;
}

void MaterialType::SetStrainMeasures(const Vector &u, const DenseMatrix& gradu, const IntegrationPoint &ip, ElementTransformation &T)
{
   int dim = u.Size();
   F.SetSize(dim);

   // deformation gradient
   F = Identity_Mat;
   F += gradu;

   // detJ
   J = F.Det();
   if (J <= 1e-15)
   {
      std::cout << "J = " << J << std::endl;
      MFEM_ABORT("Error: J <= 0");     
   }

   // E= 0.5*( gradu + graduT + graduT*gradu)
   E.SetSize(dim);
   E = gradu;
   DenseMatrix graduT(gradu);
   graduT.Transpose();
   E += graduT;
   DenseMatrix graduT_gradu(dim);
   graduT_gradu = 0.0;
   Mult(graduT,gradu,graduT_gradu);
   E += graduT_gradu;
   E *= 0.5;

   // Left Cauchy-Green Deformation Tensor b = FF^T
   // so we can output Cauchy Stress for both formulations
   b.SetSize(dim);

   Ft = F;
   Ft.Transpose();
   Mult(F,Ft,b);

   SetMaterialConstants(ip,T);

   SetAdditionalStrainMeasures(u,gradu);   

   SetAdditionalData(u,gradu);
}

void MaterialType::SetMaterialConstants(const IntegrationPoint &ip, ElementTransformation &T)
{
   MFEM_ABORT("Not implemented for base class!");
}