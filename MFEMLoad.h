#pragma once
#include <FECore/FESurfaceLoad.h>

class MFEMSolver;

class MFEMLoad : public FESurfaceLoad
{
public:
	MFEMLoad(FEModel* fem) : FESurfaceLoad(fem) {}

	virtual void Apply(MFEMSolver* solver) = 0;
};

class MFEMNeumannBC : public MFEMLoad
{
public:
	MFEMNeumannBC(FEModel* fem) : MFEMLoad(fem) {}

	void Apply(MFEMSolver* solver) override;
};
