#pragma once
#include <FECore/FESolver.h>

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

   // flag for immersed solid
   bool ImmersedSolid();
   // solid solver type
   int SolidSolverType();
   enum SolidSolverType { MFEM = 0, FEBio = 1 };
   // three field formulation
   bool ThreeFieldFormulation();

   // velocity bcs
	void AddVelDirichletBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr);
   void AddVelDirichletBC(mfem::Coefficient* coeff, mfem::Array<int>& attr, int component);
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
   void AddWindkesselBC(FEParamDouble R, FEParamDouble Rd, FEParamDouble C, FEParamDouble p0, FEParamDouble pd, mfem::Array<int>& attr);
   void AddMappedInflowBC_FromFile(std::string flowfile, bool periodic, mfem::Array<int>& attr);

   // body forces
   void AddBodyForce(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr);

   // material properties
	double FluidDensity() const;
	double KinematicViscosity() const;
	double DynamicViscosity() const;

private:
	Imp& m;

   #ifdef MFEM_USE_MPI
   void PartitionMesh();
   #endif

	DECLARE_FECORE_CLASS();
};
