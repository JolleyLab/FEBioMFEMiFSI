#pragma once
#include <FECore/FESolver.h>

class MFEMSolver : public FESolver
{
	class Imp;
public:
	MFEMSolver(FEModel* fem);
	~MFEMSolver();

	bool Init() override;

	//! Solve an analysis step
	bool SolveStep() override;

	DECLARE_FECORE_CLASS();

private:
	Imp& m;
};
