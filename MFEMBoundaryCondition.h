#pragma once
#include <FECore/FESurfaceBC.h>
#include <FECore/FEModelParam.h>

class MFEMSolver;

// Base class for boundary conditions that can be used in the MFEMSolver
class MFEMBoundaryCondition : public FESurfaceBC
{
public:
	MFEMBoundaryCondition(FEModel* fem) : FESurfaceBC(fem) {}

	// we need to override this, but I don't this this is really used here.
	void CopyFrom(FEBoundaryCondition* pbc) override { assert(false); }

	// This function must be overridden by derived classes and should be used
	// to apply the boundary condition to the MFEMSolver
	virtual void Apply(MFEMSolver* solver, int attr) = 0;
};

// Dirichlet velocity boundary condition
class MFEMVelDirichletBC : public MFEMBoundaryCondition
{
public:
	MFEMVelDirichletBC(FEModel* fem) : MFEMBoundaryCondition(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:
	FEParamDouble m_vx;
	FEParamDouble m_vy;
	FEParamDouble m_vz;

	DECLARE_FECORE_CLASS();
};

// No slip bc
class MFEMNoSlipBC : public MFEMBoundaryCondition
{
public:
	MFEMNoSlipBC(FEModel* fem) : MFEMBoundaryCondition(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:

	DECLARE_FECORE_CLASS();
};

class MFEMPrescribedVelocityComponentBC : public MFEMBoundaryCondition
{
public:
	MFEMPrescribedVelocityComponentBC(FEModel* fem) : MFEMBoundaryCondition(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:
   FEParamDouble m_vval;
   int m_dof;

	DECLARE_FECORE_CLASS();
};

// Set periodic bcs (on both velocity and pressure)
class MFEMPeriodicBC : public MFEMBoundaryCondition
{
public:
	MFEMPeriodicBC(FEModel* fem) : MFEMBoundaryCondition(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:

	DECLARE_FECORE_CLASS();
};