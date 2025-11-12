#pragma once

#include "mfem.hpp"

using namespace mfem;
using namespace std;

// classes to handle boundary conditions in updated lagrangian approach, evaluating them on initial mesh if needed

// displacement boundary condition evaluated on initial configuration
class InitConfigDispCoeff : public VectorCoefficient
{
protected:
   // spatial dimension
   int dim;

   // displacement bc defined on initial configuration
   VectorCoefficient *total_disp = nullptr; // not owned

   // original mesh node positions
   ParGridFunction *orig_nodes = nullptr;

   // Pointer to mesh
   ParMesh *pmesh = nullptr;

   // for rebuilding initial transformation
   IsoparametricTransformation *init_trans = nullptr;
   ParFiniteElementSpace *fespace = nullptr;

public:
   InitConfigDispCoeff(int vdm_) : VectorCoefficient(vdm_) {}; // default constructor

   InitConfigDispCoeff(int vdm_, VectorCoefficient *total_disp_, ParMesh *pmesh_, ParGridFunction *orig_nodes_) : VectorCoefficient(vdm_),
   total_disp(total_disp_), orig_nodes(orig_nodes_), pmesh(pmesh_)
   {
      fespace = orig_nodes->ParFESpace();
      dim = fespace->GetParMesh()->Dimension();
   };

   virtual void Eval(Vector &u, ElementTransformation &T, const IntegrationPoint &ip) override;

   // function to build mapping from parent -> initial configuration physical element
   void BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip);

   virtual void SetTime(double t_) override;

   virtual ~InitConfigDispCoeff()
   {
      if(init_trans) delete init_trans;
      if(total_disp) delete total_disp;
   }
};


// scalar initial config coefficient
// displacement boundary condition evaluated on initial configuration
class InitConfigScalarCoeff : public Coefficient
{
protected:
   // spatial dimension
   int dim;

   // Particular component of displacement bc defined on initial configuration
   Coefficient *total_disp_comp = nullptr;
   int component;

   // original mesh node positions
   ParGridFunction *orig_nodes = nullptr;

   // Pointer to mesh
   ParMesh *pmesh = nullptr;

   // for rebuilding initial transformation
   IsoparametricTransformation *init_trans = nullptr;
   ParFiniteElementSpace *fespace = nullptr;

public:
   InitConfigScalarCoeff() : Coefficient() {}; // default constructor

   // to handle setting single component
   InitConfigScalarCoeff(int component_, Coefficient *total_disp_comp_, ParMesh *pmesh_, ParGridFunction *orig_nodes_) : 
   Coefficient(), total_disp_comp(total_disp_comp_), orig_nodes(orig_nodes_), pmesh(pmesh_), component(component_)
   {
      fespace = orig_nodes->ParFESpace();
      dim = fespace->GetParMesh()->Dimension();
   };

   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) override;

   // function to build mapping from parent -> initial configuration physical element
   void BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip);

   virtual void SetTime(double t_) override;

   virtual ~InitConfigScalarCoeff()
   {
      if(init_trans) delete init_trans;
   }
};