#include "incompressible_neoHookean.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

// void IncompressibleNeoHookean::EvalPK2Stress(const IntegrationPoint &ip, ElementTransformation &T, const Vector &u, const DenseMatrix& gradu, DenseMatrix &stress)
// {
//    stress.SetSize(u.Size());
//    stress = 0.0;

//    // Compute S = mu (detC)^(-1/3) (I - 1/3 * trC Cinv) + -p*J*Cinv
   
//    // isochoric part
//    stress = Cinv;
//    stress *= 1.0/3.0 * trC;
//    stress *= -1.0;
//    stress += Identity_Mat;
//    stress *= mu->Eval(T,ip);
//    stress *= pow(detC,-1.0/3.0);

//    // volumetric part
//    DenseMatrix temp(u.Size());
//    temp = Cinv;
//    temp *= J * -1.0 * p;

//    // S = S_iso + S_vol
//    stress += temp;
// }

// void IncompressibleNeoHookean::EvalCauchyStress(const IntegrationPoint &ip, ElementTransformation &T, const Vector &u, const DenseMatrix& gradu, DenseMatrix &stress)
// {
//    stress.SetSize(u.Size());
//    stress = 0.0;

//    // Compute sigma = sigma_iso + sigma_vol
//    // = mu*J^(-5/3)*(b - 1/3 * trb I) - pI

//    // sigma_iso
//    stress = Identity_Mat;
//    stress *= b.Trace();
//    stress *= -1.0/3.0;
//    stress += b;
//    stress *= mu->Eval(T,ip);
//    stress *= pow(J,-5.0/3.0);

//    // sigma_vol
//    DenseMatrix temp(u.Size());
//    temp = Identity_Mat;
//    temp *= -1.0*p;

//    // sigma = sigma_iso + sigma_vol
//    stress += temp;
// }

// double IncompressibleNeoHookean::EvalElasticityTensor(const IntegrationPoint &ip, ElementTransformation &T, int i, int j, int k, int l)
// {
//    // coefficients
//    double mu_val = mu->Eval(T,ip);

//    if (formulation_type == form_type::TotalLagrangian)
//    {
//       double I_4th = (1.0 / 2.0)*( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) );

//       // deviatoric part
//       double C_dev = 2.0*mu_val*pow(detC,-1.0/3.0) * ( 1.0/3.0 * trC * I_4th - 1.0/3.0 * kron_delta(i,j)*Cinv(k,l) - 1.0/3.0 * kron_delta(k,l)*Cinv(i,j) + 1.0/9.0 * trC*Cinv(i,j)*Cinv(k,l));

//       // pressure part
//       double C_p = -p*(Cinv(i,j)*Cinv(k,l) - 2.0*I_4th);

//       // evaluate material elasticity tensor
//       return C_dev + C_p;
//    }
//    else if (formulation_type == form_type::UpdatedLagrangian)
//    {
//       double i_val = 0.5*(kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l) + kron_delta(j,k));

//       // deviatoric part
//       double C_dev = 2.0*mu_val*pow(J,-5.0/3.0) * ( 1.0/3.0 * b.Trace()*i_val - 1.0/3.0 * b(i,j)*kron_delta(k,l) - 1.0/3.0*kron_delta(i,j)*b(k,l) + 1.0/9.0*b.Trace()*kron_delta(i,j)*kron_delta(k,l));

//       // pressure part
//       double C_p = -p*(kron_delta(i,j)*kron_delta(k,l) - 2.0*i_val);

//       // evaluate spatial elasticity tensor
//       return C_dev + C_p;
//    }
//    else
//    {
//       MFEM_ABORT("Error: Unsupported Formulation Type");
//       return 0.0;
//    }
// }

// void IncompressibleNeoHookean::SetAdditionalStrainMeasures(const Vector &u, const DenseMatrix& gradu)
// {
//    // Compute C=F^T*F
//    DenseMatrix Ft(vdim);
//    Ft = F;
//    Ft.Transpose();
//    Mult(Ft,F,C);

//    detC = C.Det();

//    trC = C.Trace();

//    // // Compute Cinv
//    // if (Cinv_op == nullptr)
//    // {
//    //    Cinv_op = new DenseMatrixInverse(); // indicate sym pos def, 
//    // }

//    // // if (Cinv_op) delete Cinv_op;
//    // // Cinv_op = new DenseMatrixInverse(); // indicate sym pos def, 
//    // Cinv_op->SetOperator(C);
//    // Cinv_op->GetInverseMatrix(Cinv);

//    DenseMatrixInverse Cinv_op;
//    Cinv_op.SetOperator(C);
//    Cinv_op.GetInverseMatrix(Cinv);
// }

// void IncompressibleNeoHookean::SetPressure(const double &p_)
// {
//    p = p_;
// }

// double IncompressibleNeoHookean::ElasticityTensorNorm(const IntegrationPoint &ip, ElementTransformation &T)
// {
//    MFEM_ABORT("Error: TO DO");
//    return 0.0;
// }




// using namespace mfem;
// using namespace std;
// using namespace mat;

// void TestNearlyIncompressibleNeoHookean::EvalPK2Stress(const IntegrationPoint &ip, ElementTransformation &T, const Vector &u, const DenseMatrix& gradu, DenseMatrix &stress)
// {
//    stress.SetSize(u.Size());
   
//    // Compute S = mu (detC)^(-1/3) (I - 1/3 * trC Cinv) + kappa*(J-1)*J*Cinv
   
//    // isochoric part
//    stress = Cinv;
//    stress *= -1.0/3.0 * trC;
//    stress += Identity_Mat;
//    stress *= mu->Eval(T,ip) * pow(detC,-1.0/3.0);

//    // volumetric part, S = S_iso + S_vol
//    stress.Add(J * (J - 1.0) * kappa->Eval(T,ip), Cinv);
// }

// void TestNearlyIncompressibleNeoHookean::EvalCauchyStress(const IntegrationPoint &ip, ElementTransformation &T, const Vector &u, const DenseMatrix& gradu, DenseMatrix &stress)
// {
//    stress.SetSize(u.Size());

//    // Compute sigma = sigma_iso + sigma_vol
//    // = mu*J^(-5/3)*(b - 1/3 * trb I) + kappa*(J-1)I

//    // sigma_iso
//    stress = Identity_Mat;
//    stress *= -1.0/3.0 * b.Trace();
//    stress += b;
//    stress *= mu->Eval(T,ip) * pow(J,-5.0/3.0);

//    // sigma = sigma_iso + sigma_vol
//    stress.Add(kappa->Eval(T,ip) * (J - 1.0),Identity_Mat);
// }

// double TestNearlyIncompressibleNeoHookean::EvalElasticityTensor(const IntegrationPoint &ip, ElementTransformation &T, const DenseMatrix& gradu, int i, int j, int k, int l)
// {
//    // coefficients
//    double kappa_val = kappa->Eval(T,ip);
//    double mu_val = mu->Eval(T,ip);

//    if (formulation_type == form_type::TotalLagrangian)
//    {
//       double I_4th = (1.0 / 2.0)*( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) );

//       // deviatoric part
//       double C_dev = 2.0*mu_val*pow(detC,-1.0/3.0) * ( 1.0/3.0 * trC * I_4th - 1.0/3.0 * kron_delta(i,j)*Cinv(k,l) - 1.0/3.0 * kron_delta(k,l)*Cinv(i,j) + 1.0/9.0 * trC*Cinv(i,j)*Cinv(k,l));

//       // pressure part
//       double C_p = kappa_val*(J-1.0)*J*(Cinv(i,j)*Cinv(k,l) - 2.0*I_4th);

//       // volumetric part
//       double C_vol = kappa_val * J * J * (Cinv(i,j)*Cinv(k,l));

//       // evaluate material elasticity tensor
//       return C_dev + C_p + C_vol;
//    }
//    else if (formulation_type == form_type::UpdatedLagrangian)
//    {
//       double i_val = 0.5*(kron_delta(i,k)*kron_delta(j,l) + kron_delta(i,l) * kron_delta(j,k));

//       // deviatoric part
//       double C_dev = 2.0*mu_val*pow(J,-5.0/3.0) * ( 1.0/3.0 * b.Trace()*i_val - 1.0/3.0 * b(i,j)*kron_delta(k,l) - 1.0/3.0*kron_delta(i,j)*b(k,l) + 1.0/9.0*b.Trace()*kron_delta(i,j)*kron_delta(k,l));

//       // pressure part
//       double C_p = kappa_val*(J-1.0)*(kron_delta(i,j)*kron_delta(k,l) - 2.0*i_val);

//       // volumetric part
//       double C_vol = kappa_val * J *(kron_delta(i,j)*kron_delta(k,l));

//       // evaluate spatial elasticity tensor
//       return C_dev + C_p + C_vol;
//    }
//    else
//    {
//       MFEM_ABORT("Error: Unsupported Formulation Type");
//       return 0.0;
//    }
// }

// double TestNearlyIncompressibleNeoHookean::ElasticityTensorNorm(const IntegrationPoint &ip, ElementTransformation &T)
// {
//    MFEM_ASSERT(formulation_type == form_type::UpdatedLagrangian, "Error this method is only supported for UL approach");
   
//    // BELOW IS TESTING, NOT EXACTLY THE ELASTICITY TENSOR NORM!, but an approx similar to what it is for linear elasticity

//    // coefficients
//    double kappa_val = kappa->Eval(T,ip);
//    double mu_val = mu->Eval(T,ip);

//    double lambda_val = kappa_val - (2.0/3.0) * mu_val;

//    double lambda_prime = lambda_val / J;
//    double mu_prime = (mu_val - lambda_val * log(J)) / J;

//    return lambda_prime + 2.0*mu_prime;

//    // return 10.0E3*max(kappa_val,mu_val);
// }

// void TestNearlyIncompressibleNeoHookean::SetAdditionalStrainMeasures(const Vector &u, const DenseMatrix& gradu)
// {
//    // Compute C=F^T*F
//    DenseMatrix Ft(F);
//    Ft.Transpose();
//    Mult(Ft,F,C);

//    detC = C.Det();
//    trC = C.Trace();

//    DenseMatrixInverse Cinv_op;
//    Cinv_op.SetOperator(C);
//    Cinv_op.GetInverseMatrix(Cinv);
// }

