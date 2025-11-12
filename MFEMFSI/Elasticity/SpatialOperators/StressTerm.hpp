#pragma once

#include "mfem.hpp"
#include "../Materials/material_type.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

// Nonlinear integrator for stress term, (\grad(wh),\stress(d_h))
/*
   Utilizes Voight notation following Wriggers. Nonlinear FEM.
   S = {S_11 S_22 S_33 S_12 S_23 S_13}
*/
class StressIntegrator : public NonlinearFormIntegrator
{
private:
   ConstantCoefficient *coeff = nullptr;
   MaterialType *mat_type = nullptr;

   int formulation_type;

   // for building initial parent->phys element transformation
   // original mesh node positions
   ParGridFunction *orig_nodes = nullptr;
   // for rebuilding initial transformation
   IsoparametricTransformation *init_trans = nullptr;
   ParFiniteElementSpace *fespace = nullptr;

   DenseMatrix Identity_Mat;

   // temp
   Vector shape, u, vec, Ni_dx, Nk_dx, temp_vec;
   DenseMatrix dshape, dshapedx, gradu, EF, stress, gradphys_shape, F;
   double geo_val, w, Nix, Niy, Niz, Nkx, Nky, Nkz;
   DenseMatrix CBK;

   // Incompressibility
   ParGridFunction *pres = nullptr;
   bool incompressible_flag = false;

   // voigt elasticity tensor
   DenseSymmetricMatrix *D = nullptr; // not owned

public:
   StressIntegrator(MaterialType *mat_type_, int formulation_type_, ParGridFunction *orig_nodes_) : formulation_type(formulation_type_), mat_type(mat_type_), orig_nodes(orig_nodes_), coeff(nullptr)
   {
      fespace = orig_nodes->ParFESpace();

      // form identity matrix size (dim,dim);
      int dim = orig_nodes->ParFESpace()->GetParMesh()->Dimension();
      Identity_Mat.SetSize(dim);
      Identity_Mat = 0.0;
      for(int ii = 0; ii < dim; ii++)
      {
         Identity_Mat(ii,ii) = 1.0;
      }
   };

   StressIntegrator(MaterialType *mat_type_, int formulation_type_, ParGridFunction *orig_nodes_, ParGridFunction *pres_) : formulation_type(formulation_type_), mat_type(mat_type_), orig_nodes(orig_nodes_), pres(pres_), coeff(nullptr) 
   {
      fespace = orig_nodes->ParFESpace();

      // form identity matrix size (dim,dim);
      int dim = orig_nodes->ParFESpace()->GetParMesh()->Dimension();
      Identity_Mat.SetSize(dim);
      Identity_Mat = 0.0;
      for(int ii = 0; ii < dim; ii++)
      {
         Identity_Mat(ii,ii) = 1.0;
      }
   };

   StressIntegrator(MaterialType *mat_type_, int formulation_type_, ParGridFunction *orig_nodes_, ConstantCoefficient *coeff_) : formulation_type(formulation_type_), mat_type(mat_type_), orig_nodes(orig_nodes_), coeff(coeff_) 
   {
      fespace = orig_nodes->ParFESpace();

      // form identity matrix size (dim,dim);
      int dim = orig_nodes->ParFESpace()->GetParMesh()->Dimension();
      Identity_Mat.SetSize(dim);
      Identity_Mat = 0.0;
      for(int ii = 0; ii < dim; ii++)
      {
         Identity_Mat(ii,ii) = 1.0;
      }
   };

   virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &trans,
                                      const Vector &elfun,
                                      Vector &elvect);

   virtual void AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &trans,
                                     const Vector &elfun,
                                     DenseMatrix &elmat);

   enum form_type {TotalLagrangian = 0, UpdatedLagrangian = 1};

   // methods to handle building initial parent-> initial configuration physical element

   // function to build mapping from parent -> initial configuration physical element
   void BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip);

   virtual ~StressIntegrator() 
   {
      if(init_trans) delete init_trans;
   };
};