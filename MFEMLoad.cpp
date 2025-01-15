#include "MFEMLoad.h"
#include "MFEMSolver.h"
#include <math.h>

class TExact3D
{
public:
	MFEMSolver* mfem = nullptr;
	TExact3D(MFEMSolver* solver) : mfem(solver) { This = this; }

	double p_exact_3d(const mfem::Vector& xval, double t)
	{
		double x = xval(0);
		double y = xval(1);
		double z = xval(2);

		double kinvis = mfem->KinematicViscosity();
		double fluid_density = mfem->FluidDensity();

		double time_part = -fluid_density * M_PI * M_PI / 32.0 * exp(-M_PI * M_PI * t * kinvis / 2.0);

		double space_part = exp(M_PI * x / 2.0) + exp(M_PI * y / 2.0) + exp(M_PI * z / 2.0)
			+ 2.0 * exp(M_PI * (x + y) / 4.0) * cos(M_PI * (y + 2.0 * z) / 4.0) * sin(M_PI * (2.0 * x + z) / 4.0)
			+ 2.0 * exp(M_PI * (x + z) / 4.0) * cos(M_PI * (x + 2.0 * y) / 4.0) * sin(M_PI * (y + 2.0 * z) / 4.0)
			+ 2.0 * exp(M_PI * (y + z) / 4.0) * cos(M_PI * (2.0 * x + z) / 4.0) * sin(M_PI * (x + 2.0 * y) / 4.0);

		return time_part * space_part;
	}

	void traction_bc_3d(const mfem::Vector& xval, double t, mfem::Vector& tract)
	{
		double x = xval(0);
		double y = xval(1);
		double z = xval(2);

		tract.SetSize(xval.Size());
		tract = 0.0;

		double p = p_exact_3d(xval, t);

		double kinvis = mfem->KinematicViscosity();
		double dyn_viscosity = mfem->DynamicViscosity();

		tract(0) = p + (M_PI * M_PI * dyn_viscosity * exp(-(M_PI * M_PI * kinvis * t) / 4.0) * (exp((M_PI * x) / 4.0) * sin((M_PI * (y + 2.0 * z)) / 4.0) - exp((M_PI * z) / 4.0) * sin((M_PI * (x + 2.0 * y)) / 4.0))) / 8.0;
		tract(1) = 2.0 * dyn_viscosity * ((M_PI * exp(-(M_PI * M_PI * kinvis * t) / 4.0) * ((M_PI * exp((M_PI * x) / 4.0) * cos((M_PI * (y + 2.0 * z)) / 4.0)) / 4.0 + (M_PI * exp((M_PI * y) / 4.0) * cos((M_PI * (2.0 * x + z)) / 4.0)) / 2.0)) / 8.0 - (M_PI * exp(-(M_PI * M_PI * kinvis * t) / 4.0) * ((M_PI * exp((M_PI * z) / 4.0) * sin((M_PI * (x + 2.0 * y)) / 4.0)) / 2.0 - (M_PI * exp((M_PI * x) / 4.0) * cos((M_PI * (y + 2.0 * z)) / 4.0)) / 4.0)) / 8.0);
		tract(2) = 2.0 * dyn_viscosity * ((M_PI * exp(-(M_PI * M_PI * kinvis * t) / 4.0) * ((M_PI * exp((M_PI * x) / 4.0) * cos((M_PI * (y + 2.0 * z)) / 4.0)) / 2.0 + (M_PI * exp((M_PI * z) / 4.0) * cos((M_PI * (x + 2.0 * y)) / 4.0)) / 4.0)) / 8.0 - (M_PI * exp(-(M_PI * M_PI * kinvis * t) / 4.0) * ((M_PI * exp((M_PI * y) / 4.0) * sin((M_PI * (2.0 * x + z)) / 4.0)) / 2.0 - (M_PI * exp((M_PI * z) / 4.0) * cos((M_PI * (x + 2.0 * y)) / 4.0)) / 4.0)) / 8.0);
	}

	static void traction_bc_3d_static(const mfem::Vector& xval, double t, mfem::Vector& u)
	{
		This->traction_bc_3d(xval, t, u);
	}

	static TExact3D* This;
};

TExact3D* TExact3D::This = nullptr;

void MFEMNeumannBC::Apply(MFEMSolver* solver)
{
	TExact3D* tmp = new TExact3D(solver);

	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> vel_n_attr(mesh->bdr_attributes.Size());
	vel_n_attr = 0;
	vel_n_attr[4] = 1;
	solver->AddTractionBC(&TExact3D::traction_bc_3d_static, vel_n_attr);
}
