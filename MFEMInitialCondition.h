#pragma once
#include <FECore/FEInitialCondition.h>

class MFEMSolver;

namespace mfem {
	class Vector;
}

class MFEMInitialCondition : public FEInitialCondition
{
public:
	MFEMInitialCondition(FEModel* fem) : FEInitialCondition(fem) {}

	virtual void Apply(MFEMSolver* solver) = 0;
};

class MFEMInitialPressure : public MFEMInitialCondition
{
public:
	MFEMInitialPressure(FEModel* fem) : MFEMInitialCondition(fem) {}

	void Apply(MFEMSolver* solver) override;
};

class MFEMInitialVelocity : public MFEMInitialCondition
{
public:
	MFEMInitialVelocity(FEModel* fem) : MFEMInitialCondition(fem) {}

	void Apply(MFEMSolver* solver) override;
};
