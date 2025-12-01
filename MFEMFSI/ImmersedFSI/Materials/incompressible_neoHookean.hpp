#pragma once

// #include "mfem.hpp"
#include "material_type.hpp"

using namespace std;
using namespace mat;

// // incompressible neo hookean model where incompressibility is enforced via lagrange multiplier
// class IncompressibleNeoHookean : public MaterialType
// {
// protected:
//    // elastic coefficients
//    Coefficient *mu;

//    // additional strain measure
//    DenseMatrix C, Cinv;
//    double detC, trC;
//    double p;

// public:
//    IncompressibleNeoHookean(int form_type, int vdim, Coefficient &mu_) : MaterialType(form_type, vdim) 
//    {
//       mu = &mu_;

//       C.SetSize(vdim);
//       C = 0.0;

//       Cinv.SetSize(vdim);
//       Cinv = 0.0;
//    };

//    // pk2 at quad pt
//    virtual void EvalPK2Stress(const IntegrationPoint &ip, ElementTransformation &T, const Vector &u, const DenseMatrix& gradu, DenseMatrix &stress) override;

//    // cauchy stress at quad pt
//    virtual void EvalCauchyStress(const IntegrationPoint &ip, ElementTransformation &T, const Vector &u, const DenseMatrix& gradu, DenseMatrix &stress) override;

//    // method to evaluate Elasticity Tensor
//    virtual double EvalElasticityTensor(const IntegrationPoint &ip, ElementTransformation &T, int i, int j, int k, int l) override;

//    // virtual function called within SetStrainMeasures to compute additional quantities needed specific to a material model
//    virtual void SetAdditionalStrainMeasures(const Vector &u, const DenseMatrix& gradu) override;

//    // method to compute l2 norm of elasticity tensor
//    virtual double ElasticityTensorNorm(const IntegrationPoint &ip, ElementTransformation &T) override;

//    virtual void SetPressure(const double &p_) override;

//    virtual ~IncompressibleNeoHookean()
//    {
//       // delete Cinv_op;
//    };
// };

// #pragma once

// // #include "mfem.hpp"
// #include "material_type.hpp"

// using namespace mfem;
// using namespace std;
// using namespace mat;

// // nearly incompressible neo hookean model where incompressibility is enforced via an augmented strain energy
// class TestNearlyIncompressibleNeoHookean : public MaterialType
// {
// protected:
//    // elastic coefficients
//    Coefficient *kappa;
//    Coefficient *mu;

//    // additional strain measure
//    DenseMatrix C, Cinv;
//    double detC, trC;

// public:
//    TestNearlyIncompressibleNeoHookean(int form_type, int vdim, Coefficient &kappa_, Coefficient &mu_) : MaterialType(form_type, vdim) 
//    {
//       kappa = &kappa_;
//       mu = &mu_;

//       C.SetSize(vdim);
//       C = 0.0;

//       Cinv.SetSize(vdim);
//       Cinv = 0.0;
//    };

//    // pk2 at quad pt
//    virtual void EvalPK2Stress(const IntegrationPoint &ip, ElementTransformation &T, const Vector &u, const DenseMatrix& gradu, DenseMatrix &stress) override;

//    // cauchy stress at quad pt
//    virtual void EvalCauchyStress(const IntegrationPoint &ip, ElementTransformation &T, const Vector &u, const DenseMatrix& gradu, DenseMatrix &stress) override;

//    // method to evaluate Elasticity Tensor
//    virtual double EvalElasticityTensor(const IntegrationPoint &ip, ElementTransformation &T, const DenseMatrix& gradu, int i, int j, int k, int l) override;

//    // virtual function called within SetStrainMeasures to compute additional quantities needed specific to a material model
//    virtual void SetAdditionalStrainMeasures(const Vector &u, const DenseMatrix& gradu) override;

//    // method to compute l2 norm of elasticity tensor
//    virtual double ElasticityTensorNorm(const IntegrationPoint &ip, ElementTransformation &T) override;

//    TestNearlyIncompressibleNeoHookean()
//    {
//       // delete Cinv_op;
//    }
// };