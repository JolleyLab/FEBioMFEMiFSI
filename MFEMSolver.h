#pragma once
#include <FECore/FESolver.h>

namespace mfem {
	class GridFunction;
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
	mfem::GridFunction* GetCurrentPressure();
	mfem::GridFunction* GetCurrentVelocity();

	double FluidDensity() const;
	double KinematicViscosity() const;

private:
	Imp& m;

	DECLARE_FECORE_CLASS();
};
