#include "MFEMLoad.h"
#include "MFEMSolver.h"
#include "MFEMTools.h"
#include <FECore/FEMesh.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(MFEMSurfaceLoad,FESurfaceLoad)
   ADD_PARAMETER(m_backflow_stab, "backflow_stab");
   ADD_PARAMETER(m_beta, "beta");
END_FECORE_CLASS();

mfem::Mesh* MFEMSurfaceLoad::GetMesh(MFEMSolver* solver)
{
   if(solver->ImmersedSolid())
   {  
      #ifdef MFEM_USE_MPI
      return dynamic_cast<mfem::Mesh*>(solver->GetMFEMfluidmesh());
      #else
      return solver->GetMFEMfluidmesh();
      #endif
   }
   else
   {      
      #ifdef MFEM_USE_MPI
      return dynamic_cast<mfem::Mesh*>(solver->GetMFEMMesh());
      #else
      return solver->GetMFEMMesh();
      #endif
   }
}

void MFEMSurfaceLoad::CopyFrom(FEModelLoad* pml)
{
   MFEMSurfaceLoad* m_load = dynamic_cast<MFEMSurfaceLoad*>(pml);
   MFEM_ASSERT(m_load, "MFEMSurfaceLoad::CopyFrom failed");

   // FESurface* surf = &m_load->GetSurface();
   // SetSurface(surf);

   (*this).m_bdr_attr = m_load->m_bdr_attr;
   (*this).m_backflow_stab = m_load->m_backflow_stab;
   (*this).m_beta = m_load->m_beta;
}

void MFEMSurfaceLoad::Apply(MFEMSolver* solver, int attr)
{
   if (m_backflow_stab)
   {
      mfem::Mesh* mesh = GetMesh(solver);
      mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
      load_attr = 0;
      load_attr[attr] = 1;
      solver->AddBackflowStabilizationBC(m_beta, load_attr);
   }
}

BEGIN_FECORE_CLASS(MFEMNeumannBC, MFEMSurfaceLoad)
	ADD_PARAMETER(m_tx, "tx");
	ADD_PARAMETER(m_ty, "ty");
	ADD_PARAMETER(m_tz, "tz");
END_FECORE_CLASS();

void MFEMNeumannBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddTractionBC(new FEVectorCoefficient(GetFEModel(), m_tx, m_ty, m_tz), load_attr);
   MFEMSurfaceLoad::Apply(solver, attr);
}

void MFEMNeumannBC::CopyFrom(FEModelLoad* pml)
{
   MFEMSurfaceLoad::CopyFrom(pml);
   MFEMNeumannBC* m_nbc = dynamic_cast<MFEMNeumannBC*>(pml);
   MFEM_ASSERT(m_nbc, "MFEMNeumannBC::CopyFrom failed");
   CopyParameterListState(m_nbc->GetParameterList());
   (*this).m_tx = m_nbc->m_tx;
   (*this).m_ty = m_nbc->m_ty;
   (*this).m_tz = m_nbc->m_tz;
}

BEGIN_FECORE_CLASS(MFEMTractionFreeBC, MFEMSurfaceLoad)
END_FECORE_CLASS();

void MFEMTractionFreeBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddTractionFreeBC(load_attr);
   MFEMSurfaceLoad::Apply(solver, attr);
}

void MFEMTractionFreeBC::CopyFrom(FEModelLoad* pml)
{
   MFEMSurfaceLoad::CopyFrom(pml);
   MFEMTractionFreeBC* m_tbc = dynamic_cast<MFEMTractionFreeBC*>(pml);
   MFEM_ASSERT(m_tbc, "MFEMTractionFreeBC::CopyFrom failed");
   CopyParameterListState(m_tbc->GetParameterList());
}

BEGIN_FECORE_CLASS(MFEMPressureBC, MFEMSurfaceLoad)
	ADD_PARAMETER(m_p, "p")->setLongName("Pressure (Sets t = p*nor)");
END_FECORE_CLASS();

void MFEMPressureBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddPressureBC(new FEScalarCoefficient(GetFEModel(), m_p), load_attr);
   MFEMSurfaceLoad::Apply(solver, attr);
}

void MFEMPressureBC::CopyFrom(FEModelLoad* pml)
{
   MFEMSurfaceLoad::CopyFrom(pml);
   MFEMPressureBC* m_pbc = dynamic_cast<MFEMPressureBC*>(pml);
   MFEM_ASSERT(m_pbc, "MFEMPressureBC::CopyFrom failed");
   CopyParameterListState(m_pbc->GetParameterList());
   (*this).m_p = m_pbc->m_p;
}

// BEGIN_FECORE_CLASS(MFEMBackflowStabilizationBC, MFEMSurfaceLoad)
// 	ADD_PARAMETER(m_beta, "beta");
// END_FECORE_CLASS();

// MFEMBackflowStabilizationBC::MFEMBackflowStabilizationBC(FEModel* fem) : MFEMSurfaceLoad(fem)
// {
//    m_beta = 0.2;
// }

// void MFEMBackflowStabilizationBC::Apply(MFEMSolver* solver, int attr)
// {
// 	mfem::Mesh* mesh = GetMesh(solver);
// 	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
// 	load_attr = 0;
// 	load_attr[attr] = 1;
// 	solver->AddBackflowStabilizationBC(m_beta, load_attr);
// }

// void MFEMBackflowStabilizationBC::CopyFrom(FEModelLoad* pml)
// {
//    MFEMSurfaceLoad::CopyFrom(pml);
//    MFEMBackflowStabilizationBC* m_bf = dynamic_cast<MFEMBackflowStabilizationBC*>(pml);
//    MFEM_ASSERT(m_bf, "MFEMBackflowStabilizationBC::CopyFrom failed");
//    CopyParameterListState(m_bf->GetParameterList());
//    (*this).m_beta = m_bf->m_beta;
// }

BEGIN_FECORE_CLASS(MFEMResistanceBC, MFEMSurfaceLoad)
	ADD_PARAMETER(m_R, "R");
	ADD_PARAMETER(m_p0, "p_0");
END_FECORE_CLASS();

void MFEMResistanceBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddResistanceBC(m_R, m_p0, load_attr);
   MFEMSurfaceLoad::Apply(solver, attr);
}

void MFEMResistanceBC::CopyFrom(FEModelLoad* pml)
{
   MFEMSurfaceLoad::CopyFrom(pml);
   MFEMResistanceBC* m_rbc = dynamic_cast<MFEMResistanceBC*>(pml);
   MFEM_ASSERT(m_rbc, "MFEMResistanceBC::CopyFrom failed");
   CopyParameterListState(m_rbc->GetParameterList());
   (*this).m_R = m_rbc->m_R;
   (*this).m_p0 = m_rbc->m_p0;
}

BEGIN_FECORE_CLASS(MFEMWindkesselBC, MFEMSurfaceLoad)
	ADD_PARAMETER(m_R, "R");
	ADD_PARAMETER(m_Rd, "Rd");
   ADD_PARAMETER(m_C, "C");
   ADD_PARAMETER(m_p0, "p_0");
   ADD_PARAMETER(m_pd, "p_d");
END_FECORE_CLASS();

void MFEMWindkesselBC::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> load_attr(mesh->bdr_attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddWindkesselBC(m_R, m_Rd, m_C, m_p0, m_pd, load_attr);
   MFEMSurfaceLoad::Apply(solver, attr);
}

void MFEMWindkesselBC::CopyFrom(FEModelLoad* pml)
{
   MFEMSurfaceLoad::CopyFrom(pml);
   MFEMWindkesselBC* m_wbc = dynamic_cast<MFEMWindkesselBC*>(pml);
   MFEM_ASSERT(m_wbc, "MFEMWindkesselBC::CopyFrom failed");
   CopyParameterListState(m_wbc->GetParameterList());
   (*this).m_R = m_wbc->m_R;
   (*this).m_Rd = m_wbc->m_Rd;
   (*this).m_C = m_wbc->m_C;
   (*this).m_p0 = m_wbc->m_p0;
   (*this).m_pd = m_wbc->m_pd;
   (*this).m_Q0 = m_wbc->m_Q0;

}

mfem::Mesh* MFEMBodyLoad::GetMesh(MFEMSolver* solver)
{
   if(solver->ImmersedSolid())
   {  
      #ifdef MFEM_USE_MPI
      return dynamic_cast<mfem::Mesh*>(solver->GetMFEMfluidmesh());
      #else
      return solver->GetMFEMfluidmesh();
      #endif
   }
   else
   {      
      #ifdef MFEM_USE_MPI
      return dynamic_cast<mfem::Mesh*>(solver->GetMFEMMesh());
      #else
      return solver->GetMFEMMesh();
      #endif
   }
}

void MFEMBodyLoad::CopyFrom(FEModelLoad* pml)
{
   MFEMBodyLoad* m_bl = dynamic_cast<MFEMBodyLoad*>(pml);
   MFEM_ASSERT(m_bl, "MFEMBodyLoad::CopyFrom failed");
   
   (*this).m_dom_list = m_bl->m_dom_list;
}

BEGIN_FECORE_CLASS(MFEMBodyForce, MFEMBodyLoad)
	ADD_PARAMETER(m_bx, "bx");
	ADD_PARAMETER(m_by, "by");
   ADD_PARAMETER(m_bz, "bz");
END_FECORE_CLASS();

void MFEMBodyForce::Apply(MFEMSolver* solver, int attr)
{
	mfem::Mesh* mesh = GetMesh(solver);
	mfem::Array<int> load_attr(mesh->attributes.Max());
	load_attr = 0;
	load_attr[attr] = 1;
	solver->AddBodyForce(new FEVectorCoefficient(GetFEModel(), m_bx, m_by, m_bz), load_attr);
}

void MFEMBodyForce::CopyFrom(FEModelLoad* pml)
{
   MFEMBodyLoad::CopyFrom(pml);

   MFEMBodyForce* m_bf = dynamic_cast<MFEMBodyForce*>(pml);
   MFEM_ASSERT(m_bf, "MFEMBodyForce::CopyFrom failed");
   CopyParameterListState(m_bf->GetParameterList());

   (*this).m_bx = m_bf->m_bx;
   (*this).m_by = m_bf->m_by;
   (*this).m_bz = m_bf->m_bz;
}

BEGIN_FECORE_CLASS(MFEMSetNodalPressure, MFEMNodalLoad)
	ADD_PARAMETER(m_p, "p");
   ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list:pressure)");
END_FECORE_CLASS();

void MFEMSetNodalPressure::Apply(MFEMSolver* solver)
{
   FEMesh& m = GetMesh();

   feLogWarning("SetNodalPressure condition may cause issues for iFSI if node is within overlapping solid region!");

   // check to ensure pressure dof is set
   if (m_dof != 0)
   {
      feLogError("dof must be set to pressure for the SetNodalPressure condition!");
      assert(false);
   }   

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