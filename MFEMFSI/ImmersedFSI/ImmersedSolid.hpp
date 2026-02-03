#pragma once
#define IMMERSEDSOLID_VERSION 1.0

#include <map>

#include "mfem.hpp"
#include "Materials/material_type.hpp"
#include "ImmersedFSIFindPts.hpp"
#include "ImmersedSolidIntegrator.hpp"
#include "../../MFEMImSolid/MFEMNodalAverageOutputs.h"
#include "../IncompressibleFluid/ifs_bc_containers.hpp"

#ifdef MFEM_USE_MPI
#include "immersed_bc_containers.hpp"
#endif

using namespace std;
using namespace mat;

namespace immersed_solid
{

// Base class for immersed solids
class ImmersedSolid
{
public:
   // constructor
   ImmersedSolid() {};
   #ifdef MFEM_USE_MPI
   ImmersedSolid(mfem::ParMesh *solid_mesh, int order_, mfem::Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
	   mfem::ParFiniteElementSpace *vfes_, mfem::ParMesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_);
   #else
   ImmersedSolid(mfem::Mesh *solid_mesh, int order_, mfem::Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
	   mfem::FiniteElementSpace *vfes_, mfem::Mesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_);
   #endif


   // driving functions //
   virtual void Setup(double dt, double alpha_m_, double alpha_f_, double gamma_, double beta_);
   virtual void SetTimeAndDt(double time_, double dt_);
   virtual void ComputeSolidResidualAndSpread2Fluid(mfem::Vector &lambda_fsi_vec);
   virtual void Predictors(mfem::Vector &vf, mfem::Vector &vdotf);
   virtual void SetEssentialBCs(const double &time);
   virtual void SetWeakBCs(const double &time);
   virtual void Initiator();
   virtual void Correctors(mfem::Vector &del);
   virtual void Update();
   virtual void UpdatePointSearch(const mfem::GridFunction &d, const mfem::Vector &true_d);
   virtual void MoveMesh(mfem::GridFunction &total_disp_gf);
   virtual mfem::Operator &FormLHS();
   virtual void ComputeMortarIntegral()
   {
      #ifdef MFEM_USE_GSLIB
      sw_fsicmat->Start();
      finder_fluid->AssembleCouplingMatrix(true_interp_points_in_solid,vfes,dfes,dfes->GetOrdering());     
      sw_fsicmat->Stop();    
      #endif
   };

   // Get currrent gridfunctions //
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction *GetCurrentDisplacement()
   {
      return &dn_gf;
   };
   mfem::ParGridFunction *GetCurrentVelocity()
   {     
      return &ddotn_gf; 
   };
   mfem::ParGridFunction *GetCurrentAcceleration()
   {     
      return &d2dotn_gf; 
   };
   #else
   mfem::GridFunction *GetCurrentDisplacement()
   {
      return &dn_gf;
   };
   mfem::GridFunction *GetCurrentVelocity()
   {     
      return &ddotn_gf; 
   };
   mfem::GridFunction *GetCurrentAcceleration()
   {     
      return &d2dotn_gf; 
   };
   #endif   

   // Get the mesh
   #ifdef MFEM_USE_MPI
   mfem::ParMesh *GetMesh() { return pmesh; }
   #else
   mfem::Mesh *GetMesh() { return pmesh; }
   #endif

   // SETTINGS //
   // Set solution unknown ordering type (Default is byVdim ie node 1 x,y,z then node 2 x,y,z etc)
   void SetDofOrderingType(int ordering_type) { orderingtype = ordering_type; };

   // Set material type for stress
   void SetMaterialType(MaterialType *mat_type_) { mat_type = mat_type_; };

   // set formulation type
   enum form_type {TotalLagrangian = 0, UpdatedLagrangian = 1};
   void SetFormulationType(int form_)
   {
      formulation_type = form_;
   }

   void StabsScaling(mfem::GridFunction *s_gf_fluid_);

   void SetTimers(mfem::StopWatch *sw_res_,
                  mfem::StopWatch *sw_tan_,
                  mfem::StopWatch *sw_fsires_,
                  mfem::StopWatch *sw_fsitan_,
                  mfem::StopWatch *sw_fsicmat_)
   {
      sw_tan = sw_tan_;
      sw_res = sw_res_;
      sw_fsires = sw_fsires_;
      sw_fsitan = sw_fsitan_;
      sw_fsicmat = sw_fsicmat_;
   }

   // track the motion of a node in the solid mesh
   void TrackNodalDisplacement(mfem::Vector &coord, mfem::Vector &disp, double tol = 1e-12);

   #ifdef MFEM_USE_MPI
   // BCs
   /// Add a traction boundary condition
   void AddTractionBC(mfem::VectorCoefficient *coeff, mfem::Array<int> &attr);
   void AddTractionBC(VecFuncT *f, mfem::Array<int> &attr);

   /// Add a follower pressure load boundary condition, p > 0 here (for compression set load < 0)
   void AddFollowerPressureLoad(mfem::Coefficient *coeff, mfem::Array<int> &attr);
   void AddFollowerPressureLoad(ScalarFuncT *f, mfem::Array<int> &attr);
   #endif

   // restart functions
   virtual void writeRestarts(std::ofstream &sout, std::string solidmname, std::string dname, std::string ddotname, std::string d2dotname);
   virtual void readRestarts(std::ifstream &in, mfem::BaseKDTreeNodalProjection* map_d, mfem::BaseKDTreeNodalProjection* map_ddot, mfem::BaseKDTreeNodalProjection* map_d2dot, std::string solidmname, std::string dname, std::string ddotname, std::string d2dotname);
   virtual void UpdateSolidStateFromRestart();

   // Output
   // compute nodally averaged stress/strain, returned as scalar gridfunctions
   virtual mfem::GridFunction *ComputeNodallyAveragedStress(nodal_avg::StressQuantity measure, mfem::GridFunction &d_gf, bool elemental_avg = false);
   virtual mfem::GridFunction *ComputeNodallyAveragedStrain(nodal_avg::StrainQuantity measure, mfem::GridFunction &d_gf, bool elemental_avg = false);
   void BuildInitTransformation(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);

   virtual ~ImmersedSolid();

protected:
   /// Print information about the solver
   void PrintInfo();

   // default options //
   bool verbose = true;

   // flag for updated or total lagrangian framework
   int formulation_type = form_type::TotalLagrangian; // default is TL

   // // choice of ordering solution unknowns
   int orderingtype = mfem::Ordering::byNODES;

   // solid mesh   
   #ifdef MFEM_USE_MPI
   mfem::ParMesh *pmesh = nullptr; /// not owned
   mfem::ParGridFunction orig_nodes, new_nodes; // mesh nodes
   #else
   mfem::Mesh *pmesh = nullptr; /// not owned
   mfem::GridFunction orig_nodes, new_nodes; // mesh nodes
   #endif

   // solid properties
   mfem::Coefficient *init_rhos = nullptr;
   double mu_s = 0.0; // for viscoelastic solids

   // Displacement finite element collection (H^1) or L2 space for dg method
   mfem::FiniteElementCollection *dfec = nullptr;

   // Displacement (H^1)^d finite element space
   #ifdef MFEM_USE_MPI
   mfem::ParFiniteElementSpace *dfes = nullptr;
   #else
   mfem::FiniteElementSpace *dfes = nullptr;
   #endif

   // All displacement boundaries
   mfem::Array<int> disp_attr;
   // std::vector<DisplacementBC_T> disp_bcs;
   
   // essential dof list
   mfem::Array<int> disp_ess_tdof_list; // prescribed and fixed disp tdof list

   // material type
   MaterialType *mat_type = nullptr; // not owned

   // parameters
   int orderd;
   double time, dt;

   #ifdef MFEM_USE_MPI
   using NodalField = mfem::ParGridFunction;
   #else
   using NodalField = mfem::GridFunction;
   #endif

   // nodal output storage keyed by requested measure
   std::map<nodal_avg::StressQuantity, NodalField> nodal_stress_fields;
   std::map<nodal_avg::StrainQuantity, NodalField> nodal_strain_fields;

   // Operators
   mfem::SumCoefficient *density_diff = nullptr;
   #ifdef MFEM_USE_MPI
   mfem::ParNonlinearForm *iSolidNLF = nullptr;
   #else
   mfem::NonlinearForm *iSolidNLF = nullptr;
   #endif
   
   #ifdef MFEM_USE_MPI
   // Traction BCs
   mfem::ParLinearForm *Ft_form = nullptr;
   // Bookkeeping for traction bcs
   std::vector<TractionBC_T> traction_bcs;
   // All traction boundaries
   mfem::Array<int> traction_attr;
   mfem::Vector ft;

   // follower pressure load BCs
   mfem::ConstantCoefficient minusOne;
   // Bookkeeping for follower pressure loads
   std::vector<FollowerPressureLoadBC_T> fp_bcs;
   #endif

   // fluid data
   double fluid_density;
   double mu_f;
   #ifdef MFEM_USE_MPI
   mfem::ParFiniteElementSpace *vfes = nullptr; // not owned
   mfem::ParMesh *fluid_pmesh = nullptr; // not owned  
   #else
   mfem::FiniteElementSpace *vfes = nullptr; // not owned
   mfem::Mesh *fluid_pmesh = nullptr; // not owned  
   #endif    
   std::vector<AccelTerm_T> *accel_terms = nullptr;
   int my_rank = 0;

   // Solution vectors and Matrices   
   mfem::Vector res, ddotg, d2dotg, dg, ddotnp1, d2dotnp1, d2dotn, ddotn, dnp1, dn, delvdot_ib;
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction dn_gf, ddotn_gf, d2dotn_gf, curr_d2dot, curr_ddot, temp_delvdot_ib_gf, dtemp_gf, temp_del_gf;
   #else
   mfem::GridFunction dn_gf, ddotn_gf, d2dotn_gf, curr_d2dot, curr_ddot, temp_delvdot_ib_gf, dtemp_gf, temp_del_gf;
   #endif

   #ifdef MFEM_USE_MPI
   mfem::HypreParMatrix *KCt = nullptr;
   mfem::HypreParMatrix *CKCt = nullptr;
   #else
   mfem::SparseMatrix *KCt = nullptr;
   mfem::SparseMatrix *CKCt = nullptr;
   #endif   

   // time stepping
   double alpha_f, alpha_m, gamma, beta;

   // point search
   #ifdef MFEM_USE_GSLIB
   ImmersedFSIFindPts *finder_fluid = nullptr;
   #endif
   mfem::Vector interp_points_in_solid, true_interp_points_in_solid, interp_points_in_solid_tn; // fluid->solid

   // track displacement of a node
   int local_node_num = -1;
   bool find_node = true;

   //  settings
   bool total_integ = true;
   bool subtract_fluid_visc = true; // subtract artificial fluid viscous stress, on by default

   // VMS stabilization scaling
   mfem::GridFunction *s_gf_fluid = nullptr; // not owned
   bool stabs_scaling = false;
   mfem::SubMesh *solid_bdr_mesh = nullptr;
   mfem::Vector solid_bdr_nodes, solid_bdr_nodes_tn;
   mfem::Array<int> solid_bdr_tdofs;

   // timers
   mfem::StopWatch *sw_tan = nullptr;
   mfem::StopWatch *sw_res = nullptr;
   mfem::StopWatch *sw_fsires = nullptr;
   mfem::StopWatch *sw_fsitan = nullptr;
   mfem::StopWatch *sw_fsicmat = nullptr;

   // output data
   #ifdef MFEM_USE_MPI
   mfem::ParFiniteElementSpace *nodal_stress_fes = nullptr;
   #else
   mfem::FiniteElementSpace *nodal_stress_fes = nullptr;
   #endif
   mfem::FiniteElementCollection *nodal_stress_fec = nullptr;
   mfem::FiniteElementCollection *nodal_strain_fec = nullptr;
   #ifdef MFEM_USE_MPI
   mfem::ParFiniteElementSpace *nodal_strain_fes = nullptr;
   #else
   mfem::FiniteElementSpace *nodal_strain_fes = nullptr;
   #endif
   mfem::IsoparametricTransformation *init_trans = nullptr;

};

}; // namespace immersed_fsi