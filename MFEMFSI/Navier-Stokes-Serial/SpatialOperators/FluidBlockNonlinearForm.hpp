#pragma once

#include "mfem.hpp"
#include "../../Navier-Stokes/ins_bc_containers.hpp"
#include "../../Navier-Stokes/ins_solver_enum.hpp"

namespace ins_solver_serial
{

// Coefficient for tau_m = tau_SUPG = tau_PSPG //

// standard choice
// note: may need to account for nodal permutations for simplex elements in def of G dot G (see J. Liu CMAME 2020)
class tau_m_coeff : public mfem::Coefficient
{   
public:
   double kinvis;
   double CT = 4.0; // standard choice
   double CI = 36.0; // for LINEAR shape functions
   double *dt = nullptr;
   int dim;

   mfem::DenseMatrix G, K, temp;
   bool first_flag = true;

   // coefficient for weakening stab but increasing volume loss penalization
   // used in immersed fsi
   // scale = s_gf(ip) * scale_fac
   mfem::GridFunction *s_gf = nullptr;
   double scale_fac = 1.0;

   void UseConsistentTimeTerm()
   {
      consistent_time_term = true;
   };
   bool consistent_time_term = false;
   double *udot_norm = nullptr;
   double *u_norm = nullptr; // in ale case this will be norm of u - uhat

   // option for additional momentum restoring term
   void UseAdditionalMomentumRestoringTerm() { momentum_restoring_term = true; };
   void SetTol(const double &eps_) { eps = eps_; };
   bool momentum_restoring_term = false;
   double eps = std::numeric_limits<double>::epsilon();

   tau_m_coeff() : mfem::Coefficient() {};
   tau_m_coeff(double kinvis_) : mfem::Coefficient(), kinvis(kinvis_) {};
   tau_m_coeff(double kinvis_, double CI_) : mfem::Coefficient(), kinvis(kinvis_), CI(CI_) {};

   double EvalTauB(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip, const double &tau_m, const mfem::Vector &rm)
   {
      if (momentum_restoring_term)
      {
         const mfem::DenseMatrix &dxi_dx = T.InverseJacobian();
         if (first_flag)
         {
            dim = rm.Size();
            G.SetSize(dim);
            SetupMats(T);
            first_flag = false;
         }

         // apply fix for nodal permutations on triangles and tets
         if (T.mesh->GetElementType(T.ElementNo) == mfem::Element::TRIANGLE)
         {
            Mult(K,dxi_dx,temp);
            MultAtB(dxi_dx,temp,G);
         }
         else if (T.mesh->GetElementType(T.ElementNo) == mfem::Element::TETRAHEDRON)
         {
            Mult(K,dxi_dx,temp);
            MultAtB(dxi_dx,temp,G);
         }
         else
         {
            MultAtB(dxi_dx,dxi_dx,G); // G = (dxi/dx)^T (dxi/dx)
         }     

         double val = G.InnerProduct(rm,rm);
         val *= tau_m * tau_m;

         if (abs(val) < eps) val = eps; //////////////////////////////////////////
         return 1.0/sqrt(val);
      }
      else
      {
         return 0.0;
      }
   };

   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
   {
      MFEM_ABORT("not implemented");
      return 0.0;
   };

   double EvalTauM(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip, const mfem::Vector &u)
   {
      const mfem::DenseMatrix &dxi_dx = T.InverseJacobian();
      if (first_flag)
      {
         dim = u.Size();
         G.SetSize(u.Size());
         SetupMats(T);
         first_flag = false;
      }

      // apply fix for nodal permutations on triangles and tets
      if (T.mesh->GetElementType(T.ElementNo) == mfem::Element::TRIANGLE)
      {
         Mult(K,dxi_dx,temp);
         MultAtB(dxi_dx,temp,G);
      }
      else if (T.mesh->GetElementType(T.ElementNo) == mfem::Element::TETRAHEDRON)
      {
         Mult(K,dxi_dx,temp);
         MultAtB(dxi_dx,temp,G);
      }
      else
      {
         MultAtB(dxi_dx,dxi_dx,G); // G = (dxi/dx)^T (dxi/dx)
      }      

      // viscous
      double G_dot_G = G.FNorm2();      
      double val = G_dot_G * kinvis * kinvis * CI;

      // convective
      val += G.InnerProduct(u,u);

      // robustness in zero-velocity inviscid limit
      val += std::numeric_limits<double>::epsilon();

      // temporal
      if (consistent_time_term)
      {
         double omega = (*udot_norm) / (*u_norm);
         val += omega * omega;

         double temp = CT / ((*dt) * (*dt));
         double temp2 = 10.0;
      }
      else
      {
         val += CT / ((*dt) * (*dt));
      }     

      // vms scaling factor
      if (s_gf) 
      {
         MFEM_ASSERT(scale_fac >= 1.0, "VMS stabilization scaling factor must be >= 1.0!");
         double scale = scale_fac * s_gf->GetValue(T.ElementNo,ip);
         if (abs(scale) > 0.0) val *= scale;
      }

      return 1.0/sqrt(val);
   };

   mfem::DenseMatrix *GetG()
   {
      return &G;
   }

   double GetTraceG()
   {
      return G.Trace();
   }

   void SetupMats(mfem::ElementTransformation &T)
   {
      if (T.mesh->GetElementType(T.ElementNo) == mfem::Element::TRIANGLE)
      {
         K.SetSize(dim);
         K(0,0) = 2.0;
         K(0,1) = -1.0;
         K(1,0) = -1.0;
         K(1,1) = 2.0;
         K *= 1.0 / sqrt(3.0);

         temp.SetSize(dim);
      }
      else if (T.mesh->GetElementType(T.ElementNo) == mfem::Element::TETRAHEDRON)
      {
         K.SetSize(dim);
         K(0,0) = 3.0;
         K(0,1) = -1.0;
         K(0,2) = -1.0;

         K(1,0) = -1.0;
         K(1,1) = 3.0;
         K(1,2) = -1.0;

         K(2,0) = -1.0;
         K(2,1) = -1.0;
         K(2,2) = 3.0;

         K *= 1.0 / (2.0 * cbrt(2.0));

         temp.SetSize(dim);
      }
   }

   void Set_CI(int &orderv)
   {
      CI = 60.0 * pow(2.0, (double) orderv - 2.0);
      // modification following https://doi.org/10.1016/j.compfluid.2015.05.011
   }

   virtual ~tau_m_coeff() {};
};

// Coefficient tau_c = (tr G * tau_m) ^ -1 //
class nu_LSIC : public mfem::Coefficient
{
public:
   tau_m_coeff *tau_coeff = nullptr;

   nu_LSIC() : mfem::Coefficient() {};
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
   {
      MFEM_ABORT("not implemented");
      return 0.0;
   };

   virtual double Eval_Coeff(const double &tau_m)
   {
      double trG = tau_coeff->GetTraceG();

      return 1.0 / (trG * tau_m);
   };
   virtual ~nu_LSIC() {};
};

#ifdef MFEM_USE_MPI
// Special fluid nonlinear form class for FSI //
class FluidParBlockNonlinearForm : public mfem::ParBlockNonlinearForm
{
private:
   // Local to fluid mesh velocity and pressure space
   
   mfem::ParFiniteElementSpace *vfes = nullptr; // not owned
   mfem::ParFiniteElementSpace *pfes = nullptr; // not owned

   // temps
   mfem::BlockVector xs_true, ys_true, xs, ys;

   mutable mfem::BlockVector aux1_fsi;

   void MapFluidDofs2GlobalDofs(mfem::Array<int> &fluid_dofs,mfem::Array<int> &glob_dofs);
   mfem::Array<int> fluid_2_global_dofs, global_vdofs_v;

   int global_bdr_attr_size;

protected:

   /// Set of Boundary Integrators to be applied.
   mfem::Array<mfem::BlockNonlinearFormIntegrator*> boundary_integs;
   mfem::Array<mfem::Array<int>*> boundary_integs_marker; ///< Entries are not owned.


public:

   FluidParBlockNonlinearForm(mfem::Array<mfem::ParFiniteElementSpace *> &pf, mfem::Array<int> fluid_2_glob_, int glob_bdr_attr_size_) :
   mfem::ParBlockNonlinearForm(pf), vfes(pf[0]), pfes(pf[1]), fluid_2_global_dofs(fluid_2_glob_), global_bdr_attr_size(glob_bdr_attr_size_) {};

   virtual void Mult(const mfem::Vector &x, mfem::Vector &y) const;

   void MultFSI(const mfem::Vector &y, mfem::BlockVector &res);

   void GetGradientFSI(const mfem::Vector &y, mfem::Array2D<mfem::SparseMatrix *> &BlockGrad);

   const mfem::BlockVector &ProlongateBlock(const mfem::BlockVector &x) const;

   /// Adds new Boundary Integrator. Assumes ownership of @a nbfi.
   void AddBoundaryIntegratorFSI(mfem::BlockNonlinearFormIntegrator *nbfi)
   {
      boundary_integs.Append (nbfi);
      boundary_integs_marker.Append(NULL); // NULL marker means apply everywhere
   }

   /** @brief Adds new Boundary Integrator, restricted to specific boundary
       attributes.

      Assumes ownership of @a nbfi. The array @a bdr_marker is stored internally
      as a pointer to the given Array<int> object. */
   void AddBoundaryIntegratorFSI(mfem::BlockNonlinearFormIntegrator *nbfi,
                              mfem::Array<int> &bdr_marker)
   {
      boundary_integs.Append (nbfi);
      boundary_integs_marker.Append(&bdr_marker);
   }

   /// Access all the integrators added with AddBoundaryIntegrator().
   mfem::Array<mfem::BlockNonlinearFormIntegrator*> *GetBBFI() { return &boundary_integs; }
   /** @brief Access all boundary markers added with AddBoundaryIntegrator().
       If no marker was specified when the integrator was added, the
      corresponding pointer (to Array<int>) will be NULL. */
   mfem::Array<mfem::Array<int>*> *GetBBFI_Marker() { return &boundary_integs_marker; }

   virtual ~FluidParBlockNonlinearForm();
};
#endif

// INTEGRATORS //

class Fluid_Integrator : public mfem::BlockNonlinearFormIntegrator
{
private:
   // Generalized alpha
   double alpha_m = 0.0;
   double alpha_f = 0.0;
   double gamma = 0.0;
   mfem::GridFunction *vdot_gf = nullptr; // not owned
   double *dt = nullptr;

   std::vector<AccelTerm_T> *body_force_terms = nullptr;

   // temp data
   mfem::Vector shape_u, u, u_dot_gradu, shape_u_lap, u_lap;
   mfem::DenseMatrix dshape_u, dshapedx_u, gradu, EF;
   double divu;

   mfem::Vector shape_p, gradp;
   mfem::DenseMatrix dshape_p, dshapedx_p;
   double p;

   mfem::Vector Nadx, Nbdx, Nadx_sigma;
   mfem::Vector rm, rm_inertial, gradu_rm, f_vec;
   mfem::DenseMatrix fluid_stress;

   // flags for stabilization terms
   bool supg_flag = false;
   bool pspg_flag = false;
   bool vms_flag = false;
   bool lsic_flag = false;
   
   tau_m_coeff *tau_m_evaluator = nullptr;
   mfem::Coefficient *rhof = nullptr;
   mfem::Coefficient *dyn_vis = nullptr;
   nu_LSIC *mass_stab_coeff = nullptr;
   double tau_m = 0.0;
   double tau_b = 0.0;
   double tau_c = 0.0;

   bool nxx_flag = false; // flag to compute Nxx terms
   // above is only needed for p>1 to ensure proper convergence
   // but does improve accuracy for p=1 on elements other than triangles and tets

   // ALE
   mfem::VectorGridFunctionCoefficient *mesh_velo = nullptr;
   mfem::Vector u_mesh;

   // flag for viscous term form
   int viscous_term = ins_solver::ViscousTermType::StressForm;

   // computes strong form momentum residual
   void Compute_rm(mfem::Vector &rm, mfem::Vector &rm_inertial, double &mu_f, mfem::ElementTransformation &trans, const mfem::IntegrationPoint &ip, const mfem::Array<const mfem::FiniteElement *> &el);
   
public:
   Fluid_Integrator() {};

   /// Perform the local action of the NonlinearFormIntegrator
   void AssembleElementVector(const mfem::Array<const mfem::FiniteElement *> &el,
                                    mfem::ElementTransformation &trans,
                                    const mfem::Array<const mfem::Vector *> &elfun,
                                    const mfem::Array<mfem::Vector *> &elvec) override;

   /// Assemble the local gradient matrix
   void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*> &el,
                                    mfem::ElementTransformation &trans,
                                    const mfem::Array<const mfem::Vector *> &elfun,
                                    const mfem::Array2D<mfem::DenseMatrix *> &elmats) override;


   // Set relevant data
   void SetFluidProperties(mfem::Coefficient &rhof_, mfem::Coefficient *dyn_vis_)
   {
      rhof = &rhof_;
      dyn_vis = dyn_vis_;
   };

   void SetFluidData(mfem::GridFunction *vdot_gf_, std::vector<AccelTerm_T> *body_forces_, int visc_term_type = ins_solver::ViscousTermType::StressForm)
   {
      vdot_gf = vdot_gf_;
      body_force_terms = body_forces_;
      viscous_term = visc_term_type;
   };

   void SetTimeData(double &alpha_m_, double &alpha_f_, double &gamma_, double *dt_) 
   { 
      alpha_m = alpha_m_;
      alpha_f = alpha_f_;
      gamma = gamma_;
      dt = dt_; 
   };

   void SetStabilizationTerms(bool supg_flag_, bool pspg_flag_, bool vms_flag_, bool lsic_flag_, tau_m_coeff &tau_m_, nu_LSIC &m_coeff)
   {
      supg_flag = supg_flag_;
      pspg_flag = pspg_flag_;
      vms_flag = vms_flag_;
      mass_stab_coeff = &m_coeff;
      tau_m_evaluator = &tau_m_;
      lsic_flag = lsic_flag_;
   };

   void SetALEData(mfem::VectorGridFunctionCoefficient *mesh_velo_)
   {
      mesh_velo = mesh_velo_;
   };

   void NxxFlag(bool flag)
   {
      nxx_flag = flag;
   };

   virtual ~Fluid_Integrator() {};
};
};