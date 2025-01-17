#pragma once
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEModelParam.h>

class MFEMSolver;

// Base class for boundary conditions that can be used in the MFEMSolver
class MFEMBoundaryCondition : public FEBoundaryCondition
{
public:
	MFEMBoundaryCondition(FEModel* fem) : FEBoundaryCondition(fem) {}

	// we need to override this, but I don't this this is really used here.
	void CopyFrom(FEBoundaryCondition* pbc) override { assert(false); }

	// This function must be overridden by derived classes and should be used
	// to apply the boundary condition to the MFEMSolver
	virtual void Apply(MFEMSolver* solver) = 0;
};

// Dirichlet velocity boundary condition
class MFEMVelDirichletBC : public MFEMBoundaryCondition
{
public:
	MFEMVelDirichletBC(FEModel* fem) : MFEMBoundaryCondition(fem) {}

	void Apply(MFEMSolver* solver) override;

private:
	FEParamDouble m_vx;
	FEParamDouble m_vy;
	FEParamDouble m_vz;

	DECLARE_FECORE_CLASS();
};
