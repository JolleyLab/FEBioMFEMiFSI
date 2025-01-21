#pragma once
#include <FECore/FESolver.h>

#ifdef WIN32
#include "ins_solver.hpp"
#else
#include <Navier-Stokes-Serial/ins_solver.hpp>
#endif

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
	mfem::GridFunction* GetCurrentPressure();
	mfem::GridFunction* GetCurrentVelocity();

	mfem::Mesh* GetMFEMMesh();

   // velocity bcs
	void AddVelDirichletBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr);
   void AddVelDirichletBC(mfem::Coefficient* coeff, mfem::Array<int>& attr, int component);
   void AddNoSlipBC(mfem::Array<int>& attr);
   void AddNormalVelocityBC(mfem::Coefficient* coeff, mfem::Array<int>& attr);
   // to do - add mapped velocity bc

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

   // material properties
	double FluidDensity() const;
	double KinematicViscosity() const;
	double DynamicViscosity() const;

private:
	Imp& m;

	DECLARE_FECORE_CLASS();
};
