#pragma once

#include "mfem.hpp"

using namespace mfem;
using namespace std;

// Nonlinear integrator for incompressibility constraint term \int_ref q (J-1)
class IncConstraintIntegrator : public NonlinearFormIntegrator
{
private:
   int formulation_type;
   int orderd;

   ParGridFunction *u = nullptr;
   // to build disp FE
   // FiniteElement *fe_disp = nullptr;

   // for building initial parent->phys element transformation
   // original mesh node positions
   ParGridFunction *orig_nodes = nullptr;
   // for rebuilding initial transformation
   IsoparametricTransformation *init_trans = nullptr;
   ParFiniteElementSpace *fespace = nullptr;

   DenseMatrix Identity_Mat;

   // temp
   // Vector shape, u, vec;
   // DenseMatrix dshape, dshapedx, gradu, EF, stress, gradphys_shape, temp_mat, temp_mat2;

public:
   IncConstraintIntegrator(int formulation_type_, ParGridFunction *orig_nodes_, ParGridFunction *u_, int orderd_) 
   : orig_nodes(orig_nodes_), formulation_type(formulation_type_), u(u_), orderd(orderd_)
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

   enum form_type {TotalLagrangian = 0, UpdatedLagrangian = 1};

   // methods to handle building initial parent-> initial configuration physical element
   // function to build mapping from parent -> initial configuration physical element
   void BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip);

   // double ComputeJminus1(ElementTransformation &T, const IntegrationPoint &ip);

   virtual ~IncConstraintIntegrator()
   {
      // delete fe_disp;
   };
};