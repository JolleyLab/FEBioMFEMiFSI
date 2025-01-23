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
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Max());
	vel_bc_attr = 0;
	vel_bc_attr[attr] = 1;
	solver->AddVelDirichletBC(new FEVectorCoefficient(GetFEModel(), m_vx, m_vy, m_vz), vel_bc_attr);
}

BEGIN_FECORE_CLASS(MFEMNoSlipBC, MFEMBoundaryCondition)
	ADD_PARAMETER(m_vx, "vx");
	ADD_PARAMETER(m_vy, "vy");
	ADD_PARAMETER(m_vz, "vz");
END_FECORE_CLASS();

void MFEMNoSlipBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Max());
	vel_bc_attr = 0;
	vel_bc_attr[attr] = 1;
   m_val = 0.0;
   if (m_vx && m_vy && m_vz)
   {
      solver->AddNoSlipBC(vel_bc_attr);
   }
   else
   {
      if (m_vx)
      {
         solver->AddVelDirichletBC(new FEScalarCoefficient(GetFEModel(),m_val), vel_bc_attr, 1);
      }
      if (m_vy)
      {
         solver->AddVelDirichletBC(new FEScalarCoefficient(GetFEModel(),m_val), vel_bc_attr, 2);
      }
      if(m_vz)
      {
         solver->AddVelDirichletBC(new FEScalarCoefficient(GetFEModel(),m_val), vel_bc_attr, 3);
      }
   }
}

BEGIN_FECORE_CLASS(MFEMNormalVelocityBC, MFEMBoundaryCondition)
	ADD_PARAMETER(m_vn, "vn");
END_FECORE_CLASS();

void MFEMNormalVelocityBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Max());
	vel_bc_attr = 0;
	vel_bc_attr[attr] = 1;
   solver->AddNormalVelocityBC(new FEScalarCoefficient(GetFEModel(),m_vn), vel_bc_attr);
}

BEGIN_FECORE_CLASS(MFEMPeriodicBC, MFEMBoundaryCondition)
END_FECORE_CLASS();

void MFEMPeriodicBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = solver->GetMFEMMesh();
	mfem::Array<int> p_bc_attr(mesh->bdr_attributes.Max());
	p_bc_attr = 0;
	p_bc_attr[attr] = 1;
   solver->AddPeriodicBC(p_bc_attr);
}