#include "MFEMLoad.h"
#include "MFEMSolver.h"
#include "MFEMTools.h"
#include <FECore/FEMesh.h>

BEGIN_FECORE_CLASS(MFEMNeumannBC, MFEMLoad)
	ADD_PARAMETER(m_tx, "tx");
	ADD_PARAMETER(m_ty, "ty");
	ADD_PARAMETER(m_tz, "tz");
END_FECORE_CLASS();

void MFEMNeumannBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddTractionBC(new FEVectorCoefficient(GetFEModel(), m_tx, m_ty, m_tz), load_attr);
}

BEGIN_FECORE_CLASS(MFEMTractionFreeBC, MFEMLoad)
END_FECORE_CLASS();

void MFEMTractionFreeBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddTractionFreeBC(load_attr);
}

BEGIN_FECORE_CLASS(MFEMPressureBC, MFEMLoad)
	ADD_PARAMETER(m_p, "p")->setLongName("Pressure (Sets t = p*nor)");
END_FECORE_CLASS();

void MFEMPressureBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddPressureBC(new FEScalarCoefficient(GetFEModel(), m_p), load_attr);
}

BEGIN_FECORE_CLASS(MFEMBackflowStabilizationBC, MFEMLoad)
	ADD_PARAMETER(m_beta, "beta");
END_FECORE_CLASS();

MFEMBackflowStabilizationBC::MFEMBackflowStabilizationBC(FEModel* fem) : MFEMLoad(fem)
{
   m_beta = 0.2;
}

void MFEMBackflowStabilizationBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddBackflowStabilizationBC(m_beta, load_attr);
}

BEGIN_FECORE_CLASS(MFEMResistanceBC, MFEMLoad)
	ADD_PARAMETER(m_R, "R");
	ADD_PARAMETER(m_p0, "P(0)");
END_FECORE_CLASS();

void MFEMResistanceBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddResistanceBC(m_R, m_p0, load_attr);
}

BEGIN_FECORE_CLASS(MFEMWindkesselBC, MFEMLoad)
	ADD_PARAMETER(m_R, "R");
	ADD_PARAMETER(m_Rd, "Rd");
   ADD_PARAMETER(m_C, "C");
   ADD_PARAMETER(m_p0, "p_0");
   ADD_PARAMETER(m_pd, "p_d");
END_FECORE_CLASS();

void MFEMWindkesselBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddWindkesselBC(m_R, m_Rd, m_C, m_p0, m_pd, load_attr);
}

BEGIN_FECORE_CLASS(MFEMBodyForce, MFEMBodyLoad)
	ADD_PARAMETER(m_bx, "bx");
	ADD_PARAMETER(m_by, "by");
   ADD_PARAMETER(m_bz, "bz");
END_FECORE_CLASS();

void MFEMBodyForce::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> load_attr(mesh->attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddBodyForce(new FEVectorCoefficient(GetFEModel(), m_bx, m_by, m_bz), load_attr);
}

BEGIN_FECORE_CLASS(MFEMSetNodalPressure, MFEMNodalLoad)
	ADD_PARAMETER(m_p, "p");
END_FECORE_CLASS();

void MFEMSetNodalPressure::Apply(MFEMSolver* solver)
{
   FEMesh& m = GetMesh();

	// check to ensure only one node is set
   FENodeSet& nset = *m_nodeSet;
   int N = nset.Size();
   assert(N == 1);
   
   // get the nodal position in the reference state
	vec3d r = m.Node(nset[0]).m_r0;
   mfem::Vector node_pos(3);
   node_pos(0) = r.x;
   node_pos(1) = r.y;
   node_pos(2) = r.z;

   solver->SetNodalPressure(m_p, node_pos);
}