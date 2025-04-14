#include "MFEMInitialCondition.h"
#include "MFEMSolver.h"
#ifdef WIN32
#include <mfem/fem/gridfunc.hpp>
#else
#include "mfem.hpp"
#endif
#include "MFEMTools.h"

BEGIN_FECORE_CLASS(MFEMInitialPressure, MFEMInitialCondition)
	ADD_PARAMETER(m_pressure, "pressure");
END_FECORE_CLASS();

void MFEMInitialPressure::Apply(MFEMSolver* solver)
{
	mfem::GridFunction* p_ic = solver->GetCurrentPressure();
	FEScalarCoefficient p_coeff(GetFEModel(), m_pressure);
	p_coeff.SetTime(0.0);
	p_ic->ProjectCoefficient(p_coeff);
}

BEGIN_FECORE_CLASS(MFEMInitialVelocity, MFEMInitialCondition)
	ADD_PARAMETER(m_vx, "vx");
	ADD_PARAMETER(m_vy, "vy");
	ADD_PARAMETER(m_vz, "vz");
END_FECORE_CLASS();

void MFEMInitialVelocity::Apply(MFEMSolver* solver)
{
	mfem::GridFunction* u_ic = solver->GetCurrentVelocity();
	FEVectorCoefficient u_coeff(GetFEModel(), m_vx, m_vy, m_vz);
	u_coeff.SetTime(0.0);
	u_ic->ProjectCoefficient(u_coeff);
}
