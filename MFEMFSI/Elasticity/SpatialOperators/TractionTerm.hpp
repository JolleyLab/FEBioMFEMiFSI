#pragma once

#include "mfem.hpp"

using namespace mfem;
using namespace std;

/** Class for boundary integration of L(v) := (g, v), where
    f=(f1,...,fn) and v=(v1,...,vn). */
class VectorBoundaryLFIntegrator_Ref : public LinearFormIntegrator
{
private:
   Vector shape, vec;
   VectorCoefficient &Q;

   // // for building initial parent->phys element transformation
   // // original mesh node positions
   // ParGridFunction *orig_nodes = nullptr;
   // // for rebuilding initial transformation
   // IsoparametricTransformation init_trans;
   // ParFiniteElementSpace *fespace = nullptr;

   // for building initial parent->phys element transformation
   // original mesh node positions
   ParGridFunction *orig_nodes = nullptr;
   // for rebuilding initial transformation
   IsoparametricTransformation *init_trans = nullptr;
   ParFiniteElementSpace *fespace = nullptr;

   int formulation_type;

public:
   /// Constructs a boundary integrator with a given VectorCoefficient QG
   VectorBoundaryLFIntegrator_Ref(VectorCoefficient &QG, int formulation_type_, ParGridFunction *orig_nodes_) : Q(QG), formulation_type(formulation_type_), orig_nodes(orig_nodes_) 
   {
      fespace = orig_nodes->ParFESpace();

      // init_trans.ElementNo = -100000; // random value
   };

   /** Given a particular boundary Finite Element and a transformation (Tr)
       computes the element boundary vector, elvect. */
   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       Vector &elvect);

   // For DG spaces
   virtual void AssembleRHSElementVect(const FiniteElement &el,
                                       FaceElementTransformations &Tr,
                                       Vector &elvect);

   using LinearFormIntegrator::AssembleRHSElementVect;

   // function to build mapping from parent -> initial configuration physical element
   void BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip, const FiniteElement *curr_elem);

   enum form_type {TotalLagrangian = 0, UpdatedLagrangian = 1};

   virtual ~VectorBoundaryLFIntegrator_Ref() {};
};