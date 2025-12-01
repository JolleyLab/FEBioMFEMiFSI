#pragma once

#include "mfem.hpp"
#include "../ImmersedFSI/Materials/material_type.hpp"
#include "../Navier-Stokes/ins_bc_containers.hpp"

using namespace mat;

namespace immersed_solid
{

// INTEGRATORS //

// Nonlinear integrator for immersed solid equations
/*
   Utilizes Voight notation following Wriggers. Nonlinear FEM. for internal work term
   S = {S_11 S_22 S_33 S_12 S_23 S_13}
*/
class Immersed_Solid_Integrator : public mfem::NonlinearFormIntegrator
{
private:
	mfem::ConstantCoefficient *coeff = nullptr;
   MaterialType *mat_type = nullptr;
   mfem::Coefficient *rhos0_coeff = nullptr;
   double rhof;
   int formulation_type;

   // for building initial parent->phys element transformation
   // original mesh node positions
   mfem::GridFunction *orig_nodes = nullptr;
   // for rebuilding initial transformation
   mfem::IsoparametricTransformation *init_trans = nullptr; // owned
   mfem::FiniteElementSpace *fespace = nullptr;

   mfem::DenseMatrix Identity_Mat;

   // time stepping parameters
   double beta;
   double gamma;
   double *dt = nullptr;
   mfem::GridFunction *a_gf = nullptr;
   mfem::GridFunction *v_gf = nullptr;
   mfem::Vector inertial_force;

   // gen alpha
   double alpha_m;
   double alpha_f;
   bool gen_alpha = false;

   // body forces
   std::vector<AccelTerm_T> *bf_terms = nullptr;
   mfem::Vector bf_vec;

   // temp
   mfem::Vector shape, u, vec, Ni_dx, Nk_dx, temp_vec;
   mfem::DenseMatrix dshape, dshapedx, gradu, EF, stress, gradphys_shape, F;
   double geo_val, w, Nix, Niy, Niz, Nkx, Nky, Nkz;
   mfem::DenseMatrix CBK;

   // artificial fluid viscous stress
   mfem::Vector v_nodes, e_dNdx, dNidx_j, dNdx_i, dNdx_k,FinvT_dNdX;
   mfem::Array<int> v_vdofs;
   mfem::DenseMatrix VF,dNd_currx, e_v,Finv,FinvT;
   double mu_f;

   // voigt elasticity tensor
   mfem::DenseSymmetricMatrix *D = nullptr; // not owned

   // bool to determine whether to subtract off fluid stress
   bool subtract_fluid_visc_stress = true;
   
public:

   Immersed_Solid_Integrator(MaterialType *mat_type_, int formulation_type_, mfem::GridFunction *orig_nodes_, mfem::GridFunction *curr_a_, mfem::GridFunction *curr_v_, double *dt_, double beta_, double gamma_, double alpha_m_, double alpha_f_, std::vector<AccelTerm_T> *bf_terms_, mfem::Coefficient *init_rhos_, double &rhof_, double &mu_f_)
   : formulation_type(formulation_type_), mat_type(mat_type_), orig_nodes(orig_nodes_), coeff(nullptr), a_gf(curr_a_), v_gf(curr_v_), dt(dt_), beta(beta_), alpha_m(alpha_m_), alpha_f(alpha_f_), bf_terms(bf_terms_), rhos0_coeff(init_rhos_), rhof(rhof_), mu_f(mu_f_), gamma(gamma_)
   {
      fespace = orig_nodes->FESpace();

      // form identity matrix size (dim,dim);
      int dim = orig_nodes->FESpace()->GetMesh()->Dimension();
      Identity_Mat.SetSize(dim);
      Identity_Mat = 0.0;
      for(int ii = 0; ii < dim; ii++)
      {
         Identity_Mat(ii,ii) = 1.0;
      }
   };

   Immersed_Solid_Integrator(MaterialType *mat_type_, int formulation_type_, mfem::GridFunction *orig_nodes_, mfem::GridFunction *curr_a_, mfem::GridFunction *curr_v_, double *dt_, double beta_, double gamma_, double alpha_m_, double alpha_f_, std::vector<AccelTerm_T> *bf_terms_, mfem::Coefficient *init_rhos_, double &rhof_, double &mu_f_, bool visc_)
   : formulation_type(formulation_type_), mat_type(mat_type_), orig_nodes(orig_nodes_), coeff(nullptr), a_gf(curr_a_), v_gf(curr_v_), dt(dt_), beta(beta_), alpha_m(alpha_m_), alpha_f(alpha_f_), bf_terms(bf_terms_), rhos0_coeff(init_rhos_), rhof(rhof_), mu_f(mu_f_), gamma(gamma_), subtract_fluid_visc_stress(visc_)
   {
      fespace = orig_nodes->FESpace();

      // form identity matrix size (dim,dim);
      int dim = orig_nodes->FESpace()->GetMesh()->Dimension();
      Identity_Mat.SetSize(dim);
      Identity_Mat = 0.0;
      for(int ii = 0; ii < dim; ii++)
      {
         Identity_Mat(ii,ii) = 1.0;
      }
   };

   /// Perform the local action of the NonlinearFormIntegrator
   virtual void AssembleElementVector(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, const mfem::Vector &elfun, mfem::Vector &elvect) override;

   /// Assemble the local gradient matrix
   virtual void AssembleElementGrad(const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, const mfem::Vector &elfun, mfem::DenseMatrix &elmat) override;
   void AssembleMaterialStiffness(int &si,int &sk, const int &dim, const int &nd, double &af, mfem::DenseMatrix &elmat);

   enum form_type {TotalLagrangian = 0, UpdatedLagrangian = 1};

   // function to build mapping from parent -> initial configuration physical element
   void BuildInitTransformation(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);

   virtual ~Immersed_Solid_Integrator() 
   {
      if(init_trans) delete init_trans;
   };
};

} // namespace immersed_solid