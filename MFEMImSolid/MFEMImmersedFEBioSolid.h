#pragma once
#include "../MFEMFSI/ImmersedFSI/ImmersedSolid.hpp"
#include "../MFEMUtilities/MFEMTools.h"
#include "../MFEMLinAlg/MFEMGlobalTypes.h"
#include "../MFEMCore/MFEMLoad.h"
#include "mfem.hpp"
#include <vector>


using namespace immersed_solid;

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
   void SetupFEDomains();
   void Predictors(mfem::Vector &vf, mfem::Vector &vdotf) override;
   void SetWeakBCs(const double &time) override {}; // handled by FEBio
   void Initiator() override;
   void Correctors(mfem::Vector &del) override;
   void Update() override;
   void UpdateMaterialPointDataForPostProcessing(); // material point data to tn

   // residual
   void ComputeSolidResidualAndSpread2Fluid(mfem::Vector &lambda_fsi_vec) override;

   // tangent
   void ReformStiffness();
   mfem::Operator &FormLHS() override;

   // set data
   void SetSubMeshMaps(mfem::TransferMap *submeshsolid_to_global_d_, mfem::TransferMap* global_to_submeshsolid_d_, mfem::GridFunction *ddot_global_gf_);

   // restart functions
   void UpdateSolidStateFromRestart() override;
   void UpdateFEBioDataFromRestart();

   // output
   mfem::GridFunction *ComputeNodallyAveragedStress(nodal_avg::StressQuantity measure, mfem::GridFunction &d_gf, bool elemental_avg = false) override;
   mfem::GridFunction *ComputeNodallyAveragedStrain(nodal_avg::StrainQuantity measure, mfem::GridFunction &d_gf, bool elemental_avg = false) override;

   // Derived measures computed from already computed nodal component fields
   mfem::GridFunction *ComputeDerivedNodalStressMeasureFromComponents(nodal_avg::StressQuantity measure, bool elemental_avg = false);
   mfem::GridFunction *ComputeDerivedNodalStrainMeasureFromComponents(nodal_avg::StrainQuantity measure, bool elemental_avg = false);

   // void BuildGaussPointSpace();
   // void UpdateGaussPointFields();
   // mfem::QuadratureFunction *GetGaussPointStressComponent(nodal_avg::StressQuantity comp);
   // mfem::QuadratureFunction *GetGaussPointStrainComponent(nodal_avg::StrainQuantity comp);

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

   // Recovered nodal fields
   #ifdef MFEM_USE_MPI
   using NodalField = mfem::ParGridFunction;
   #else
   using NodalField = mfem::GridFunction;
   #endif
   // Nodal projection-based recovered fields
   std::map<nodal_avg::StressQuantity, NodalField> nodalproj_stress_fields;
   std::map<nodal_avg::StrainQuantity, NodalField> nodalproj_strain_fields;
   // Elemental-averaged recovered fields
   std::map<nodal_avg::StressQuantity, NodalField> elemental_stress_fields;
   std::map<nodal_avg::StrainQuantity, NodalField> elemental_strain_fields;

   // // Gauss-point output storage
   // mfem::QuadratureSpace *solid_qspace = nullptr; // owned
   // int solid_qorder = 0;
   // std::map<nodal_avg::StressQuantity, mfem::QuadratureFunction*> gp_stress_fields; // owned
   // std::map<nodal_avg::StrainQuantity, mfem::QuadratureFunction*> gp_strain_fields; // owned

   int PickQuadratureOrderForSolidMesh();
};


