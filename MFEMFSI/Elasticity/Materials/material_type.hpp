#pragma once

#include "mfem.hpp"

using namespace std;

namespace mat
{

// Base class for different constitutive models
class MaterialType
{
protected:
   // flag for UL vs TL
   int formulation_type;
   int vdim;

   // identity matrix
   mfem::DenseMatrix Identity_Mat;
   mfem::DenseMatrix temp_stress;

   // Strain measures
   mfem::DenseMatrix F, Ft; // deformation gradient
   double J = 0; // detF
   mfem::DenseMatrix b; // Left Cauchy-Green Deformation Tensor
   mfem::DenseMatrix E; // Green-Lagrange Strain Tensor

   // Voigt notation elasticity tensor
   // vdim == 3
   mfem::DenseSymmetricMatrix C_voigt; // Sym 6x6 matrix
   mfem::Vector temp_vec;
   double d;

   // vdim == 2 (plane strain)
   mfem::Vector temp_vec_2;
   mfem::DenseMatrix C2d_voight;

   // vector of material constants
   mfem::Vector mat_consts;
   // ordered like in constitutive equation c0, c1, c2, etc

public:

   MaterialType() : formulation_type(-1), vdim(0) {}; // default constructor
   MaterialType(int form_type, int vdim_) : formulation_type(form_type), vdim(vdim_) 
   {
      // form identity matrix size (dim,dim);
      Identity_Mat.Diag(1.0,vdim);

      if (vdim == 3)
      {
         temp_vec.SetSize(21); // vector to store data for C_voigt (sym 6x6 matrix)
      }

      temp_stress.SetSize(vdim);
      temp_stress = 0.0;

      d = (double) vdim;
   };

   void SetParameters(int form_type, int vdim_)
   {
      formulation_type = form_type;
      vdim = vdim_;
   };

   double GetJ() 
   { 
      MFEM_VERIFY(abs(J) >= 1.0E-16, "Error: J = 0, or SetStrainMeasures must be called prior to this!");
      return J; 
   };

   mfem::DenseMatrix &GetF()
   {
      return F;
   }

   // method to evaluate PK2 Stress
   virtual void EvalPK2Stress(mfem::DenseMatrix &stress);

   // method to evaluate PK1 Stress given PK2
   virtual void EvalPK1Stress(mfem::DenseMatrix &stress);

   // method to evaluate Cauchy Stress given PK1 or directly evaluate Cauchy stress
   virtual void EvalCauchyStress(mfem::DenseMatrix &stress);

   // method to evaluate Elasticity Tensor
   virtual double EvalElasticityTensor(int i, int j, int k, int l);

   mfem::DenseSymmetricMatrix * GetElasticityTensorVoigt()
   {
      C_voigt.ClearExternalData();

      // form the symmetric matrix
      // row 0
      temp_vec(0) = (*this).EvalElasticityTensor(0,0,0,0);
      temp_vec(1) = (*this).EvalElasticityTensor(0,0,1,1);
      temp_vec(2) = (*this).EvalElasticityTensor(0,0,2,2);
      temp_vec(3) = (*this).EvalElasticityTensor(0,0,0,1);
      temp_vec(4) = (*this).EvalElasticityTensor(0,0,1,2);
      temp_vec(5) = (*this).EvalElasticityTensor(0,0,0,2);

      // row 1
      temp_vec(6) = (*this).EvalElasticityTensor(1,1,1,1);
      temp_vec(7) = (*this).EvalElasticityTensor(1,1,2,2);
      temp_vec(8) = (*this).EvalElasticityTensor(1,1,0,1);
      temp_vec(9) = (*this).EvalElasticityTensor(1,1,1,2);
      temp_vec(10) = (*this).EvalElasticityTensor(1,1,0,2);

      // row 2
      temp_vec(11) = (*this).EvalElasticityTensor(2,2,2,2);
      temp_vec(12) = (*this).EvalElasticityTensor(2,2,0,1);
      temp_vec(13) = (*this).EvalElasticityTensor(2,2,1,2);
      temp_vec(14) = (*this).EvalElasticityTensor(2,2,0,2);

      // row 3
      temp_vec(15) = (*this).EvalElasticityTensor(0,1,0,1);
      temp_vec(16) = (*this).EvalElasticityTensor(1,2,0,1);
      temp_vec(17) = (*this).EvalElasticityTensor(0,2,0,1);

      // row 4
      temp_vec(18) = (*this).EvalElasticityTensor(1,2,1,2);
      temp_vec(19) = (*this).EvalElasticityTensor(1,2,0,2);

      // row 5
      temp_vec(20) = (*this).EvalElasticityTensor(0,2,0,2);
      
      C_voigt.UseExternalData(temp_vec.GetData(),6); // second arg is n for nxn

      return &C_voigt;
   }

   void Get2DElasticityTensorVoigt(mfem::DenseMatrix &C)
   {
      MFEM_VERIFY(vdim == 2, "Error: this function is only supported in 2D");
      C.SetSize(vdim + 1); // 3x3 matrix
      C(0,0) = (*this).EvalElasticityTensor(0,0,0,0); // 1111
      C(0,1) = (*this).EvalElasticityTensor(0,0,1,1); // 1122
      C(0,2) = (*this).EvalElasticityTensor(0,0,0,1); // 1112

      C(1,0) = (*this).EvalElasticityTensor(1,1,0,0); // 2211
      C(1,1) = (*this).EvalElasticityTensor(1,1,1,1); // 2222
      C(1,2) = (*this).EvalElasticityTensor(1,1,0,1); // 2212

      C(2,0) = (*this).EvalElasticityTensor(0,1,0,0); // 1211
      C(2,1) = (*this).EvalElasticityTensor(0,1,1,1); // 1222
      C(2,2) = (*this).EvalElasticityTensor(0,1,0,1); // 1212
   }

   // function to evaluate Kronecker Delta
   double kron_delta(int i, int j);

   // function to update strain measures at the current integration point
   void SetStrainMeasures(const mfem::Vector &u, const mfem::DenseMatrix& gradu, const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T);

   // virtual function called within SetStrainMeasures to compute additional quantities needed specific to a material model
   virtual void SetAdditionalStrainMeasures(const mfem::Vector &u, const mfem::DenseMatrix& gradu) {};

   // virtual function called within SetStrainMeasures to compute additional data
   virtual void SetAdditionalData(const mfem::Vector &u, const mfem::DenseMatrix& gradu) {};

   // virtual function to setup material constants array and compute material constants
   virtual void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T);

   // additional material constant functions for immersed FSI approach
   mfem::Vector &GetMaterialConstantsArray() { return mat_consts; };
   void SetMaterialConstantsArray(mfem::Vector &mat_)
   { 
      mat_consts = mat_;
      ResetMaterialConstants();
   };
   virtual void ResetMaterialConstants() { MFEM_ABORT("Not implemented for base class"); };
   virtual double GetKappa() { MFEM_ABORT("Not implemented for base class"); };

   // for incompressible materials
   virtual void SetPressure(const double &p_) {};

   // method to compute l2 norm of elasticity tensor
   virtual double ElasticityTensorNorm();
   
   enum form_type {TotalLagrangian = 0, UpdatedLagrangian = 1};

   virtual ~MaterialType() {};
};

}; // namespace mat