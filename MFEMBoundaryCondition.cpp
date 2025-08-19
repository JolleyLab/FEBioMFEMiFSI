#include "MFEMBoundaryCondition.h"
#include "MFEMSolver.h"
#include "MFEMTools.h"

mfem::Mesh* MFEMBoundaryCondition::GetMesh(MFEMSolver* solver)
{
   if(solver->ImmersedSolid())
   {  
      #ifdef iFSI_USE_MPI
      return dynamic_cast<mfem::Mesh*>(solver->GetMFEMfluidmesh());
      #else
      return solver->GetMFEMfluidmesh();
      #endif
   }
   else
   {      
      #ifdef iFSI_USE_MPI
      return dynamic_cast<mfem::Mesh*>(solver->GetMFEMMesh());
      #else
      return solver->GetMFEMMesh();
      #endif
   }
}

void MFEMBoundaryCondition::CopyFrom(FEBoundaryCondition* pbc)
{
   MFEMBoundaryCondition* m_bc = dynamic_cast<MFEMBoundaryCondition*>(pbc);
   MFEM_ASSERT(m_bc, "MFEMBoundaryCondition::CopyFrom failed");
   
   // FESurface* surf = m_bc->GetSurface();
   // SetSurface(surf);

   (*this).m_bdr_attr = m_bc->m_bdr_attr;
}

BEGIN_FECORE_CLASS(MFEMVelDirichletBC, MFEMBoundaryCondition)
	ADD_PARAMETER(m_vx, "vx");
	ADD_PARAMETER(m_vy, "vy");
	ADD_PARAMETER(m_vz, "vz");
END_FECORE_CLASS();

void MFEMVelDirichletBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Max());
	vel_bc_attr = 0;
	vel_bc_attr[attr] = 1;
	solver->AddVelDirichletBC(new FEVectorCoefficient(GetFEModel(), m_vx, m_vy, m_vz), vel_bc_attr);
}

void MFEMVelDirichletBC::CopyFrom(FEBoundaryCondition* pbc)
{
   MFEMBoundaryCondition::CopyFrom(pbc);
   MFEMVelDirichletBC* m_vbc = dynamic_cast<MFEMVelDirichletBC*>(pbc);
   MFEM_ASSERT(m_vbc, "MFEMVelDirichletBC::CopyFrom failed");
   CopyParameterListState(m_vbc->GetParameterList());
   (*this).m_vx = m_vbc->m_vx;
   (*this).m_vy = m_vbc->m_vy;
   (*this).m_vz = m_vbc->m_vz;
}

BEGIN_FECORE_CLASS(MFEMNoSlipBC, MFEMBoundaryCondition)
END_FECORE_CLASS();

void MFEMNoSlipBC::Apply(MFEMSolver* solver, int attr)
{
   mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Max());
	vel_bc_attr = 0;
	vel_bc_attr[attr] = 1;
   solver->AddNoSlipBC(vel_bc_attr);
}

void MFEMNoSlipBC::CopyFrom(FEBoundaryCondition* pbc)
{
   MFEMBoundaryCondition::CopyFrom(pbc);
   MFEMNoSlipBC* m_noslip = dynamic_cast<MFEMNoSlipBC*>(pbc);
   MFEM_ASSERT(m_noslip, "MFEMNoSlipBC::CopyFrom failed");
   CopyParameterListState(m_noslip->GetParameterList());
}

BEGIN_FECORE_CLASS(MFEMPrescribedVelocityComponentBC, MFEMBoundaryCondition)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list:velocity)");
   ADD_PARAMETER(m_vval, "value");
END_FECORE_CLASS();

void MFEMPrescribedVelocityComponentBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Max());
	vel_bc_attr = 0;
	vel_bc_attr[attr] = 1;
   solver->AddVelDirichletBC(new FEScalarCoefficient(GetFEModel(),m_vval), vel_bc_attr, m_dof);
}

void MFEMPrescribedVelocityComponentBC::CopyFrom(FEBoundaryCondition* pbc)
{
   MFEMBoundaryCondition::CopyFrom(pbc);
   MFEMPrescribedVelocityComponentBC* m_pvbc = dynamic_cast<MFEMPrescribedVelocityComponentBC*>(pbc);
   MFEM_ASSERT(m_pvbc, "MFEMPrescribedVelocityComponentBC::CopyFrom failed");
   CopyParameterListState(m_pvbc->GetParameterList());
   (*this).m_dof = m_pvbc->m_dof;
   (*this).m_vval = m_pvbc->m_vval;
}

BEGIN_FECORE_CLASS(MFEMPeriodicBC, MFEMBoundaryCondition)
END_FECORE_CLASS();

void MFEMPeriodicBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> p_bc_attr(mesh->bdr_attributes.Max());
	p_bc_attr = 0;
	p_bc_attr[attr] = 1;
   solver->AddPeriodicBC(p_bc_attr);
}

void MFEMPeriodicBC::CopyFrom(FEBoundaryCondition* pbc)
{
   MFEMBoundaryCondition::CopyFrom(pbc);
   MFEMPeriodicBC* m_periodic = dynamic_cast<MFEMPeriodicBC*>(pbc);
   MFEM_ASSERT(m_periodic, "MFEMPeriodicBC::CopyFrom failed");
   CopyParameterListState(m_periodic->GetParameterList());
}

BEGIN_FECORE_CLASS(MFEMMappedInflowBC, MFEMBoundaryCondition)
	ADD_PARAMETER(m_flowfile, "file_name");
	ADD_PARAMETER(m_periodic, "periodic");
END_FECORE_CLASS();

void MFEMMappedInflowBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> vel_bc_attr(mesh->bdr_attributes.Max());
	vel_bc_attr = 0;
	vel_bc_attr[attr] = 1;
   solver->AddMappedInflowBC_FromFile(m_flowfile,m_periodic,vel_bc_attr);
}

void MFEMMappedInflowBC::CopyFrom(FEBoundaryCondition* pbc)
{
   MFEMBoundaryCondition::CopyFrom(pbc);
   MFEMMappedInflowBC* m_vbc = dynamic_cast<MFEMMappedInflowBC*>(pbc);
   MFEM_ASSERT(m_vbc, "MFEMMappedInflowBC::CopyFrom failed");
   CopyParameterListState(m_vbc->GetParameterList());
   (*this).m_flowfile = m_vbc->m_flowfile;
   (*this).m_periodic = m_vbc->m_periodic;
}