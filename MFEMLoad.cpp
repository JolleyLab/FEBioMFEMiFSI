#include "MFEMLoad.h"
#include "MFEMSolver.h"
#include "MFEMTools.h"

BEGIN_FECORE_CLASS(MFEMNeumannBC, MFEMLoad)
	ADD_PARAMETER(m_tx, "tx");
	ADD_PARAMETER(m_ty, "ty");
	ADD_PARAMETER(m_tz, "tz");
END_FECORE_CLASS();

void MFEMNeumannBC::Apply(MFEMSolver* solver)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> vel_n_attr(mesh->bdr_attributes.Size());
	vel_n_attr = 0;
	vel_n_attr[4] = 1;
	solver->AddTractionBC(new FEVectorCoefficient(GetFEModel(), m_tx, m_ty, m_tz), vel_n_attr);
}
