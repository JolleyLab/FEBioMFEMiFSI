#pragma once

#include "mfem.hpp"
#include "../Materials/material_type.hpp"
#include "../bc_containers.hpp"

using namespace mfem;
using namespace std;
using namespace mat;

namespace dyn_elast
{
   
// class for nonlinear residual operator for newton solver
class ResidualOpElasticity: public Operator
{
private:
   HypreParMatrix *Mmat = nullptr; // not owned
   ParNonlinearForm *Stress_form = nullptr; // not owned
   Array<int> ess_tdofs;

   // owned
   mutable HypreParMatrix *Jacobian = nullptr;
   mutable HypreParMatrix *stress_mat = nullptr;
   mutable Vector temp; // size d

   bool dg_method = false;

   ParGridFunction *curr_v_gf = nullptr;
   bool robin_bc_flag = false;

   bool static_problem = false;

public:
   // constructor
   ResidualOpElasticity(HypreParMatrix *Mmat_, ParNonlinearForm *Stress_form_, Array<int> ess_tdofs_, bool static_problem_) :
   Operator(Mmat_->Height()) , Mmat(Mmat_), Stress_form(Stress_form_), ess_tdofs(ess_tdofs_), static_problem(static_problem_) 
   {
      temp.SetSize(Mmat_->Height());
      temp = 0.0;
   };

   // Compute the residual
   virtual void Mult(const Vector &k, Vector &y) const; // needed since virtual
   void Mult_PC(const Vector &d, const Vector &v, const Vector &a, const Vector &ft, const Vector &fb, Vector &res) const;

   // Compute the gradient
   virtual Operator &GetGradient(const Vector &k) const; // needed since virtual
   Operator &GetGradient_PC(const Vector &d, const Vector &v, const Vector &a, const double dt, const double am, const double af) const;

   // set dg method flag
   void Use_DG_Method()
   {
      MFEM_ABORT("TO DO");
      dg_method = true;
   };

   void RobinBCs(ParGridFunction &gf)
   {
      robin_bc_flag = true;
      curr_v_gf = &gf;
   }

   virtual ~ResidualOpElasticity();
};

}