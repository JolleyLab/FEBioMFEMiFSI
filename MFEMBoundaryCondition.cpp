#include "MFEMBoundaryCondition.h"
#include "MFEMSolver.h"
#include "MFEMTools.h"

BEGIN_FECORE_CLASS(MFEMVelDirichletBC, MFEMBoundaryCondition)
	ADD_PARAMETER(m_vx, "vx");
	ADD_PARAMETER(m_vy, "vy");
	ADD_PARAMETER(m_vz, "vz");
END_FECORE_CLASS();

void MFEMVelDirichletBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Size());
	vel_bc_attr = 0;
	vel_bc_attr[attr] = 1;
	solver->AddVelDirichletBC(new FEVectorCoefficient(GetFEModel(), m_vx, m_vy, m_vz), vel_bc_attr);
}
