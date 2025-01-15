#include "MFEMBoundaryCondition.h"
#include "MFEMSolver.h"
#include "mfem.hpp"
#include <math.h>

class VelExact3D
{
public:
	static VelExact3D* This;

	MFEMSolver* mfem;
	VelExact3D(MFEMSolver* solver) : mfem(solver) { This = this; }

	void vel_exact_3d(const mfem::Vector& xval, double t, mfem::Vector& u)
	{
		double x = xval(0);
		double y = xval(1);
		double z = xval(2);

		double kinvis = mfem->KinematicViscosity();

		double c = M_PI * exp(-M_PI * M_PI * kinvis * t / 4.0);

		u(0) = -1.0 * c / 4.0 * (exp(M_PI * x / 4.0) * sin(M_PI * (y + 2.0 * z) / 4.0) + exp(M_PI * z / 4.0) * cos(M_PI * (x + 2.0 * y) / 4.0));
		u(1) = -1.0 * c / 4.0 * (exp(M_PI * y / 4.0) * sin(M_PI * (2.0 * x + z) / 4.0) + exp(M_PI * x / 4.0) * cos(M_PI * (y + 2.0 * z) / 4.0));
		u(2) = -1.0 * c / 4.0 * (exp(M_PI * z / 4.0) * sin(M_PI * (x + 2.0 * y) / 4.0) + exp(M_PI * y / 4.0) * cos(M_PI * (2.0 * x + z) / 4.0));
	}

	static void vel_exact_3d_static(const mfem::Vector& xval, double t, mfem::Vector& u)
	{
		This->vel_exact_3d(xval, t, u);
	}
};

VelExact3D* VelExact3D::This = nullptr;

void MFEMVelDirichletBC::Apply(MFEMSolver* solver)
{
	VelExact3D* tmp = new VelExact3D(solver);

	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Size());
	vel_bc_attr = 1;
	vel_bc_attr[4] = 0;
	solver->AddVelDirichletBC(&VelExact3D::vel_exact_3d_static, vel_bc_attr);
}
