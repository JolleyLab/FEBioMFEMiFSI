#include "MFEMInitialCondition.h"
#include "MFEMSolver.h"
#include <mfem.hpp>
#include <mfem/fem/gridfunc.hpp>
#include <math.h>

class PTemp {
public:
	MFEMSolver* solver = nullptr;
	PTemp(MFEMSolver* psolver) : solver(psolver) {}

	double operator () (const mfem::Vector& xval, double t)
	{
		double x = xval(0);
		double y = xval(1);
		double z = xval(2);

		double fluid_density = solver->FluidDensity();
		double kinvis = solver->KinematicViscosity();

		double time_part = -fluid_density * M_PI * M_PI / 32.0 * exp(-M_PI * M_PI * t * kinvis / 2.0);

		double space_part = exp(M_PI * x / 2.0) + exp(M_PI * y / 2.0) + exp(M_PI * z / 2.0)
			+ 2.0 * exp(M_PI * (x + y) / 4.0) * cos(M_PI * (y + 2.0 * z) / 4.0) * sin(M_PI * (2.0 * x + z) / 4.0)
			+ 2.0 * exp(M_PI * (x + z) / 4.0) * cos(M_PI * (x + 2.0 * y) / 4.0) * sin(M_PI * (y + 2.0 * z) / 4.0)
			+ 2.0 * exp(M_PI * (y + z) / 4.0) * cos(M_PI * (2.0 * x + z) / 4.0) * sin(M_PI * (x + 2.0 * y) / 4.0);

		return time_part * space_part;
	}
};

void MFEMInitialPressure::Apply(MFEMSolver* solver)
{
	mfem::GridFunction* p_ic = solver->GetCurrentPressure();
	std::function<double(const mfem::Vector&, double)> TDF = PTemp(solver);
	mfem::FunctionCoefficient p_excoeff(TDF);
	p_excoeff.SetTime(0.0);
	p_ic->ProjectCoefficient(p_excoeff);
}

class VTemp {
public:
	MFEMSolver* solver = nullptr;
	VTemp(MFEMSolver* psolver) : solver(psolver) {}

	void operator () (const mfem::Vector& xval, double t, mfem::Vector& u)
	{
		double x = xval(0);
		double y = xval(1);
		double z = xval(2);

		double kinvis = solver->KinematicViscosity();

		double c = M_PI * exp(-M_PI * M_PI * kinvis * t / 4.0);

		u(0) = -1.0 * c / 4.0 * (exp(M_PI * x / 4.0) * sin(M_PI * (y + 2.0 * z) / 4.0) + exp(M_PI * z / 4.0) * cos(M_PI * (x + 2.0 * y) / 4.0));
		u(1) = -1.0 * c / 4.0 * (exp(M_PI * y / 4.0) * sin(M_PI * (2.0 * x + z) / 4.0) + exp(M_PI * x / 4.0) * cos(M_PI * (y + 2.0 * z) / 4.0));
		u(2) = -1.0 * c / 4.0 * (exp(M_PI * z / 4.0) * sin(M_PI * (x + 2.0 * y) / 4.0) + exp(M_PI * y / 4.0) * cos(M_PI * (2.0 * x + z) / 4.0));
	}
};

void MFEMInitialVelocity::Apply(MFEMSolver* solver)
{
	mfem::GridFunction* u_ic = solver->GetCurrentVelocity();
	std::function<void(const mfem::Vector&, double, mfem::Vector&)> TDF = VTemp(solver);
	mfem::VectorFunctionCoefficient u_excoeff(3, TDF);
	u_excoeff.SetTime(0.0);
	u_ic->ProjectCoefficient(u_excoeff);
}
