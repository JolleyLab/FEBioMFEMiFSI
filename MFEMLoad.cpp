#include "MFEMLoad.h"
#include "MFEMSolver.h"
#include "MFEMTools.h"

BEGIN_FECORE_CLASS(MFEMNeumannBC, MFEMLoad)
	ADD_PARAMETER(m_tx, "tx");
	ADD_PARAMETER(m_ty, "ty");
	ADD_PARAMETER(m_tz, "tz");
END_FECORE_CLASS();

void MFEMNeumannBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> load_attr(mesh->bdr_attributes.Size());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddTractionBC(new FEVectorCoefficient(GetFEModel(), m_tx, m_ty, m_tz), load_attr);
}
