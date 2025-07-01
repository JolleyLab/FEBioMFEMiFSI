#pragma once
#include <ImmersedFSI-Serial/immersed_solid.hpp>
#include "MFEMTools.h"
#include "MFEMGlobalTypes.h"
#include "MFEMLoad.h"
#include "mfem.hpp"
#include <vector>


using namespace immersed_solid_serial;

class ImmersedFEBioSolid : public ImmersedSolid
{
public:
   // constructors
   ImmersedFEBioSolid() : ImmersedSolid() {};
   #ifdef MFEM_USE_MPI
   ImmersedFEBioSolid(FEModel* fem_, mfem::ParMesh *solid_mesh, int order_, mfem::Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
      mfem::ParFiniteElementSpace *vfes_, mfem::ParMesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_);
   #else
   ImmersedFEBioSolid(FEModel* fem_, mfem::Mesh *solid_mesh, int order_, mfem::Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
      mfem::FiniteElementSpace *vfes_, mfem::Mesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_);
   #endif

   // driving functions
   void Setup(double dt, double alpha_m_, double alpha_f_, double gamma_, double beta_) override;
   void Predictors(mfem::Vector &vf, mfem::Vector &vdotf) override;
   void SetWeakBCs(const double &time) override {}; // handled by FEBio
   void Initiator() override;
   void Correctors(mfem::Vector &del) override;
   void Update() override;

   // residual
   void ComputeSolidResidualAndSpread2Fluid(mfem::Vector &lambda_fsi_vec) override;

   // tangent
   void ReformStiffness();
   mfem::Operator &FormLHS() override;

   // set data
   void SetSubMeshMaps(mfem::TransferMap *submeshsolid_to_global_d_, mfem::TransferMap* global_to_submeshsolid_d_, mfem::GridFunction *ddot_global_gf_);

   // destructor
   virtual ~ImmersedFEBioSolid();

public:
   // FEBio model
   FEModel *fem = nullptr;
   
   // Global residual vector
   MFEMGlobalVector *m_R = nullptr; // owned
   std::vector<double> R;
   std::vector<double> Fr;

   // Global matrix
   MFEMGlobalMatrix *m_K = nullptr; // owned

   // Model loads
   // surface loads
   bool immersed_surface_loads = false;   

   // body loads
   std::vector<MFEMBodyForce*> m_body_forces;

   // Submesh Transfer Maps
   mfem::TransferMap *submeshsolid_to_global_d = nullptr;
   mfem::TransferMap *global_to_submeshsolid_d = nullptr;
   mfem::GridFunction *ddot_global_gf = nullptr;
   mfem::GridFunction *d2dot_global_gf = nullptr; // owned
   mfem::GridFunction *dnp1_global_gf = nullptr; // owned
   mfem::GridFunction delvdot_global_gf;
   
   // Data
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction ddotnp1_gf, d2dotnp1_gf, dnp1_gf;
   #else
   mfem::GridFunction ddotnp1_gf, d2dotnp1_gf, dnp1_gf;
   #endif

   #ifdef MFEM_USE_MPI
   mfem::OperatorHandle pK;
   #endif

private:
   // flag for recomputing sparse matrix profile
   bool reshape_stiffness = true; 
};


