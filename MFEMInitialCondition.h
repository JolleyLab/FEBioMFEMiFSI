#pragma once
#include <FECore/FEInitialCondition.h>

class MFEMSolver;

// Base class for initial conditions that can be used in the MFEMSolver
class MFEMInitialCondition : public FEInitialCondition
{
public:
	MFEMInitialCondition(FEModel* fem) : FEInitialCondition(fem) {}

	virtual void Apply(MFEMSolver* solver) = 0;
};

// Class that applies an initial pressure
class MFEMInitialPressure : public MFEMInitialCondition
{
public:
	MFEMInitialPressure(FEModel* fem) : MFEMInitialCondition(fem) {}

	void Apply(MFEMSolver* solver) override;

private:
	FEParamDouble m_pressure;

	DECLARE_FECORE_CLASS();
};

// Class that applies an initial velocity
class MFEMInitialVelocity : public MFEMInitialCondition
{
public:
	MFEMInitialVelocity(FEModel* fem) : MFEMInitialCondition(fem) {}

	void Apply(MFEMSolver* solver) override;

private:
	FEParamDouble m_vx;
	FEParamDouble m_vy;
	FEParamDouble m_vz;

	DECLARE_FECORE_CLASS();
};
