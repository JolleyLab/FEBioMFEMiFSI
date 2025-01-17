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

	void AddVelDirichletBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr);

	void AddTractionBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr);

	double FluidDensity() const;
	double KinematicViscosity() const;
	double DynamicViscosity() const;

private:
	Imp& m;

	DECLARE_FECORE_CLASS();
};
