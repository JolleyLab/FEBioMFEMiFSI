#include "MFEMInitialCondition.h"
#include "MFEMSolver.h"
#include <mfem.hpp>
#include <mfem/fem/gridfunc.hpp>
#include <FECore/FEMaterialPoint.h>

// This class allows a FEParamDouble to be used in a mfem::FunctionCoefficient
class PTemp {
public:
	FEParamDouble& m_param;
	PTemp(FEParamDouble& param) : m_param(param) {}

	double operator () (const mfem::Vector& xval)
	{
		double x = xval(0);
		double y = xval(1);
		double z = xval(2);
		FEMaterialPoint pt;
		pt.m_r0 = vec3d(x, y, z);
		return m_param(pt);
	}
};

BEGIN_FECORE_CLASS(MFEMInitialPressure, MFEMInitialCondition)
	ADD_PARAMETER(m_pressure, "pressure");
END_FECORE_CLASS();

void MFEMInitialPressure::Apply(MFEMSolver* solver)
{
	mfem::GridFunction* p_ic = solver->GetCurrentPressure();
	std::function<double(const mfem::Vector&)> TDF = PTemp(m_pressure);
	mfem::FunctionCoefficient p_excoeff(TDF);
	p_excoeff.SetTime(0.0);
	p_ic->ProjectCoefficient(p_excoeff);
}

// This class allows a set of 3 FEParamDouble to be used in a mfem::VectorFunctionCoefficient
class VTemp {
public:
	FEParamDouble& m_vx;
	FEParamDouble& m_vy;
	FEParamDouble& m_vz;
	VTemp(FEParamDouble& vx, FEParamDouble& vy, FEParamDouble& vz, MFEMSolver* psolver) : m_vx(vx), m_vy(vy), m_vz(vz) {}

	void operator () (const mfem::Vector& xval, mfem::Vector& u)
	{
		double x = xval(0);
		double y = xval(1);
		double z = xval(2);
		FEMaterialPoint pt;
		pt.m_r0 = vec3d(x, y, z);
		u(0) = m_vx(pt);
		u(1) = m_vy(pt);
		u(2) = m_vz(pt);
	}
};

BEGIN_FECORE_CLASS(MFEMInitialVelocity, MFEMInitialCondition)
	ADD_PARAMETER(m_vx, "vx");
	ADD_PARAMETER(m_vy, "vy");
	ADD_PARAMETER(m_vz, "vz");
END_FECORE_CLASS();

void MFEMInitialVelocity::Apply(MFEMSolver* solver)
{
	mfem::GridFunction* u_ic = solver->GetCurrentVelocity();
	std::function<void(const mfem::Vector&, mfem::Vector&)> TDF = VTemp(m_vx, m_vy, m_vz, solver);
	mfem::VectorFunctionCoefficient u_excoeff(3, TDF);
	u_excoeff.SetTime(0.0);
	u_ic->ProjectCoefficient(u_excoeff);
}
