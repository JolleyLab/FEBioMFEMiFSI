#pragma once
#include <FECore/FESolver.h>

#include "../MFEMUtilities/MFEMParImport.h"
#include "../MFEMUtilities/MFEMParExport.h"

namespace mfem
{  
   class GridFunction;
   class Mesh;
   class VectorCoefficient;
   template <typename T> class Array;
   class Coefficient;
   class Vector;
   class VectorCoefficient;
   #ifdef MFEM_USE_MPI
   class ParGridFunction;
   class ParMesh;
   #endif
}

namespace ifs_solver
{
   class IncompressibleFluidSolver;
}

class MFEMSolver : public FESolver
{
	class Imp;
public:
	MFEMSolver(FEModel* fem);
	~MFEMSolver();

	bool Init() override;

	//! Solve an analysis step
	bool SolveStep() override;

public:
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction* GetCurrentPressure();
	mfem::ParGridFunction* GetCurrentVelocity();

	mfem::ParMesh* GetMFEMMesh();
   mfem::ParMesh* GetMFEMfluidmesh();
   mfem::ParMesh* GetMFEMsolidmesh();
   #else
   mfem::GridFunction* GetCurrentPressure();
	mfem::GridFunction* GetCurrentVelocity();

	mfem::Mesh* GetMFEMMesh();
   mfem::Mesh* GetMFEMfluidmesh();
   mfem::Mesh* GetMFEMsolidmesh();
   #endif

   ifs_solver::IncompressibleFluidSolver* GetFlowSolver();

   // flag for immersed solid
   bool ImmersedSolid();
   // solid solver type
   int SolidSolverType();
   enum SolidSolverType { MFEM = 0, FEBio = 1 };
   // three field formulation
   bool ThreeFieldFormulation();

   // velocity bcs
	void AddVelDirichletBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr, bool zero_bdr = false);
   void AddVelDirichletBC(mfem::Coefficient* coeff, mfem::Array<int>& attr, int component, bool zero_bdr = false);
   void AddNoSlipBC(mfem::Array<int>& attr);
   void AddPrescribedVelocityComponentBC(mfem::Coefficient* coeff, mfem::Array<int>& attr, int dof);
   // to do - add mapped velocity bc

   // Pressure conditions
   void SetNodalPressure(FEParamDouble pval, mfem::Vector &node_pos);

   // traction bcs
	void AddTractionBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr);
   void AddTractionFreeBC(mfem::Array<int>& attr);
   void AddPressureBC(mfem::Coefficient* coeff, mfem::Array<int>& attr);

   // periodic bcs
   void AddPeriodicBC(mfem::Array<int>& attr);

   // cardiovascular bcs
   void AddBackflowStabilizationBC(FEParamDouble beta, mfem::Array<int>& attr);
   void AddResistanceBC(FEParamDouble R, FEParamDouble p0, mfem::Array<int>& attr);
   void AddWindkesselBC(FEParamDouble R, FEParamDouble Rd, FEParamDouble C, FEParamDouble p0, FEParamDouble pd, mfem::Array<int>& attr, bool compute_initial_pressure = false);
   void AddMappedInflowBC_FromFile(std::string flowfile, bool periodic, int profile, mfem::Array<int>& attr);

   // body forces
   void AddBodyForce(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr);

   // material properties
	double FluidDensity() const;
	double KinematicViscosity() const;
	double DynamicViscosity() const;

   // time stepping
   void AdjustTimeStep();

   // Print timing data
   void PostSolve();

private:
	Imp& m;

   // Workflow functions
   bool CreateMesh();
   bool ExtractFluidAndSolidMeshes();
   bool SetupFluid();
   void SetFluidConditions(bool apply);
   bool SetupSolid();
   void SetSolidConditions(bool apply);
   void SetupSolvers();
   void PostProcessOutput();
   void FinalizeSolid();
   void RestartFiles();
   void UpdateSolidOutputs();

   #ifdef MFEM_USE_MPI
   void PartitionMesh();
   void ClearSerialMeshData();

   // helper function for pre-partitioning of meshes
   void SavePartitionedMesh(const std::string &prefix, mfem::Mesh &mesh, int np,
                         int *partitioning = nullptr);
   mfem::ParMesh* LoadParMesh(const std::string &prefix);
   #endif

   // XML writing and reading for parallel meshes //
   MFEMParImport par_import;
   MFEMParExport par_export;   

	DECLARE_FECORE_CLASS();
};
