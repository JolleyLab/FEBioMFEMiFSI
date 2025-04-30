#include "MFEMSolver.h"

#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>

#include "MFEMImmersedFEBioSolid.h"
#include "MFEMMesh.h"
#include "MFEMInitialCondition.h"
#include "MFEMBoundaryCondition.h"
#include "MFEMLoad.h"
#include "MFEMMaterial.h"
#include "FEBioLinearSolver.h"
#include "ParFEModel.h"

#include <Navier-Stokes-Serial/ins_solver.hpp>
#include <Navier-Stokes/ins_solver_enum.hpp>
#include <ImmersedFSI-Serial/immersed_solid.hpp>
#include <Elasticity/Materials/materials.hpp>

#include <typeinfo>

using namespace ins_solver_serial;
using namespace std;
using namespace std::placeholders;
using namespace immersed_solid_serial;
using namespace mat;
using namespace mfem;

#ifndef MFEM_USE_MPI
// we override this so we can use the FEBio linear solvers
class FEBio_INS_Solver : public INS_Solver
{
public:
	FEBio_INS_Solver(
		FEModel* fem,
		mfem::Mesh* mesh,
		int orderv,
		int orderp,
		double kin_vis,
		double fluid_density_,
		int ordering_type_ = mfem::Ordering::byNODES,
		bool immersed_flag = false) : INS_Solver(mesh, orderv, orderp, kin_vis, fluid_density_, ordering_type_, immersed_flag) 
	{
		m_fem = fem;
	}

	void AllocDirectSolver() override
	{
		solver = new FEBioLinearSolver(m_fem);
	}

private:
	FEModel* m_fem;
};
#endif

class MFEMSolver::Imp
{
public:
	static Imp* This;
	Imp() { This = this; }

public:
	MFEMMesh* mesh = nullptr; // MFEM mesh
   // transfer maps from submesh to global mfem mesh
   TransferMap* submeshfluid_to_global_v = nullptr; 
   TransferMap* submeshfluid_to_global_p = nullptr; 
   TransferMap* submeshsolid_to_global_d = nullptr;
   TransferMap* global_to_submeshsolid_d = nullptr;
   // global spaces and gridfunctions
   FiniteElementSpace* vfespace_global = nullptr;
   FiniteElementSpace* pfespace_global = nullptr;
   GridFunction u_gf_global, p_gf_global, d_gf_global, vib_gf_global;

   // fluid
   SubMesh* fluid_mesh = nullptr; // fluid mesh
	INS_Solver* flowsolver = nullptr; // MFEM flow solver
   #ifdef MFEM_USE_MPI
   ParGridFunction* p_gf = nullptr;
   ParGridFunction* u_gf = nullptr;
   #else
	GridFunction* p_gf = nullptr; // velocity grid function
	GridFunction* u_gf = nullptr; // pressure grid function
   #endif
	ParaViewDataCollection* pvdc = nullptr; // paraview data collection
   #ifdef MFEM_USE_MPI
   ParMesh* fluid_pmesh = nullptr;
   #endif

   // immersed solid
   SubMesh* solid_mesh = nullptr; // immersed solid mesh
   #ifdef MFEM_USE_MPI
   ParGridFunction *d_gf = nullptr; // immersed solid displacement grid function
   ParGridFunction *vib_gf = nullptr; // immersed solid velocity grid function
   #else
   GridFunction *d_gf = nullptr; // immersed solid displacement grid function
   GridFunction *vib_gf = nullptr; // immersed solid velocity grid function
   #endif
   ParaViewDataCollection *sdc = nullptr; // immersed solid paraview data collection
   int solidsolvertype = 0; // 0 = MFEM, 1 = FEBio (default is MFEM)
   ImmersedFEBioSolid *immersed_febio_solid = nullptr; // not owned
   bool three_field = false; // flag for 3-field formulation
   #ifdef MFEM_USE_MPI
   ParMesh* solid_pmesh = nullptr;
   #endif

public:
   // Material properties
   double fluid_density;
   double kinvis;
   double dyn_viscosity;

	// Discretization options
   double orderv = 1.0;
   double orderp = 1.0;
   bool stabilized_LSIC = false; // turned on automatically using VMS
   bool stabilized_VMS = true;
   int stab_type = ins_solver::StabilizationType::VMS;
   int orderingtype= Ordering::byVDIM;

   int nonlinear_form = ins_solver::NonlinearFormType::Convective;
   int viscous_term_type = ins_solver::ViscousTermType::StressForm;
   bool LaplaceForm = false;
   bool total_integ = true;
   bool zero_mean_pressure = false;

   // Time Stepping Options
   double rhoinf = 0.5;
   int step_save = 1; // frequency of storing data

   // Nonlinear Solver Options
   int minIt = 3;
   int maxIt = 10;
   double atol = 1.0e-6;
   double rtol = 1.0e-3;

   // Linear Solver Options
   // int lin_minIt;
   int lin_maxIt = 200;
   double lin_atol = 1.0e-6;
   double lin_rtol = 1.0e-3;
   bool direct_solver = false;
   // PARDISO = OpenMP (serial code only), CPARDISO = MPI (parallel code only)
   enum LinearSolverType { UMFPACK = 0, PARDISO = 1, FGMRES = 2, CPARDISO = 3 }; 
   int linsolvertype = LinearSolverType::UMFPACK;

   // Preconditioner Options
   int prec_type = ins_solver::PreconditionerType::SIMPLE;
   enum PrecInnerSolverType {Iterative = 0, Direct = 1}; // direct using UMFPACK
   int prec_inner_solver_type = PrecInnerSolverType::Iterative;
   bool direct_solver_prec = false;
   double prec_Ainv_atol = 1.0e-6;
   double prec_Ainv_rtol = 1.0e-3;
   double prec_Sinv_atol = 1.0e-6;
   double prec_Sinv_rtol = 1.0e-3;

   // immersed solid
   int orderd = orderv;
   MaterialType *mat_type = nullptr;
   bool febio_material = false;
   
   // int form_type = ImmersedSolid::form_type::TotalLagrangian; // either work
   int form_type = ImmersedSolid::form_type::UpdatedLagrangian; // either work, but FEBio computes virtual interal work integrals in spatial frame

   bool subtract_viscous_stress = true;
   double stabs_scale = 1.0;
   bool immersed_solid = false;
   int vol_type = NearlyIncompressibleMaterial::volumetric_energy::ST91;
   // TO DO make vol_type a drop down option with nearly incompressible materials

public:
   // Custom stream buffer for writing to the log file
   CustomStreamBuffer customBuffer;

public:
   #ifdef MFEM_USE_MPI
   ParFEModel* par_fem = nullptr;
   #endif
};

MFEMSolver::Imp* MFEMSolver::Imp::This = nullptr;

// TODO: Add solver parameters
BEGIN_FECORE_CLASS(MFEMSolver, FESolver)
   // Discretization Options
   BEGIN_PARAM_GROUP("Discretization");
      ADD_PARAMETER(m.orderv, FE_RANGE_GREATER_OR_EQUAL(1), "velocity_order");
      ADD_PARAMETER(m.orderp, FE_RANGE_GREATER_OR_EQUAL(1), "pressure_order");
      ADD_PARAMETER(m.stabilized_VMS, "VMS");
      // drop down menu for stabilization type, for now either VMS or unstabilized
      // ADD_PARAMETER(m.stabilized_LSIC, "LSIC");
      ADD_PARAMETER(m.LaplaceForm, "Laplace_form");
      ADD_PARAMETER(m.zero_mean_pressure, "Zero_mean_pressure");
   END_PARAM_GROUP();

   // Immersed Solid Options
   BEGIN_PARAM_GROUP("Immersed Solid Options");
      ADD_PARAMETER(m.immersed_solid, "Immersed_solid");
      ADD_PARAMETER(m.solidsolvertype, "Solid_solver_type")->setEnums("MFEM\0FEBio\0");
      ADD_PARAMETER(m.subtract_viscous_stress, "Subtract_viscous_stress");
      ADD_PARAMETER(m.three_field, "Three_field_formulation");
      ADD_PARAMETER(m.stabs_scale, FE_RANGE_NOT_EQUAL(0.0), "VMS_stabilization_scaling");
   END_PARAM_GROUP();

   // Time Stepping Options
   BEGIN_PARAM_GROUP("Time Stepping");
      ADD_PARAMETER(m.rhoinf, FE_RANGE_CLOSED(0.0,1.0), "rho_inf");
      ADD_PARAMETER(m.step_save, FE_RANGE_GREATER_OR_EQUAL(1), "step_save");
   END_PARAM_GROUP();

   // // Nonlinear Solver Options
   BEGIN_PARAM_GROUP("Nonlinear Solver");
      ADD_PARAMETER(m.atol, FE_RANGE_GREATER_OR_EQUAL(0.0), "atol");
      ADD_PARAMETER(m.rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol");
      ADD_PARAMETER(m.minIt, FE_RANGE_GREATER_OR_EQUAL(1), "minIt");
      ADD_PARAMETER(m.maxIt, FE_RANGE_GREATER_OR_EQUAL(1), "maxIt");
   END_PARAM_GROUP();

   // Linear Solver Options
   BEGIN_PARAM_GROUP("Linear Solver");
      ADD_PARAMETER(m.linsolvertype, "solver_type")->setEnums("UMFPack\0PARDISO\0FGMRES\0CPARDISO\0");
      ADD_PARAMETER(m.lin_atol, FE_RANGE_GREATER_OR_EQUAL(0.0), "lin_atol");
      ADD_PARAMETER(m.lin_rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "lin_rtol");
      // ADD_PARAMETER(m.lin_minIt, FE_RANGE_GREATER_OR_EQUAL(1), "lin_minIt");
      ADD_PARAMETER(m.lin_maxIt, FE_RANGE_GREATER_OR_EQUAL(1), "lin_maxIt");
   END_PARAM_GROUP();

   // Preconditioner Options
   BEGIN_PARAM_GROUP("Preconditioner");
      // set preconditioner drop down menu, for now SIMPLE
      ADD_PARAMETER(m.prec_inner_solver_type, "inner_solver_type")->setEnums("Iterative\0Direct\0");
      ADD_PARAMETER(m.prec_Ainv_atol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ainv_atol");
      ADD_PARAMETER(m.prec_Ainv_rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ainv_rtol");
      ADD_PARAMETER(m.prec_Sinv_atol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Sinv_atol");
      ADD_PARAMETER(m.prec_Sinv_rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Sinv_rtol");
   END_PARAM_GROUP();
END_FECORE_CLASS();

MFEMSolver::MFEMSolver(FEModel* fem) : FESolver(fem), m(*(new MFEMSolver::Imp))
{
   // default values created in Imp class

   #ifdef MFEM_USE_MPI
   // block log output not on root process
   if (!Mpi::Root()) fem->BlockLog();
   // mfem::out << "\n------------------------\n";
   // mfem::out << "Rank = " << mfem::Mpi::WorldRank() << std::endl;
   // mfem::out << "------------------------\n";
   #endif
}

MFEMSolver::~MFEMSolver()
{
   // Restore the original buffer
   mfem::out.rdbuf(m.customBuffer.originalBuffer);

   delete m.submeshfluid_to_global_v;
   delete m.submeshfluid_to_global_p;
   delete m.submeshsolid_to_global_d;
   delete m.global_to_submeshsolid_d;
   delete m.vfespace_global;
   delete m.pfespace_global;

   #ifndef MFEM_USE_MPI
	delete m.mesh;
   delete m.fluid_mesh;
   #endif
	delete m.flowsolver;
	delete m.pvdc;

   #ifndef MFEM_USE_MPI
   delete m.solid_mesh;
   #endif
   if(m.febio_material == false) delete m.mat_type;
   delete m.sdc;

   #ifdef MFEM_USE_MPI
   delete m.fluid_pmesh;
   delete m.solid_pmesh;
   delete m.par_fem;
   #endif
}

#ifdef MFEM_USE_MPI
mfem::ParGridFunction* MFEMSolver::GetCurrentPressure()
{
	return m.flowsolver->GetCurrentPressure();
}

mfem::ParGridFunction* MFEMSolver::GetCurrentVelocity()
{
	return m.flowsolver->GetCurrentVelocity();
}

mfem::ParMesh* MFEMSolver::GetMFEMMesh()
{
	return m.fluid_pmesh;
}

mfem::ParMesh* MFEMSolver::GetMFEMfluidmesh()
{
	return m.fluid_pmesh;
}

mfem::ParMesh* MFEMSolver::GetMFEMsolidmesh()
{
	return m.solid_pmesh;
}

#else
mfem::GridFunction* MFEMSolver::GetCurrentPressure()
{
	return m.flowsolver->GetCurrentPressure();
}

mfem::GridFunction* MFEMSolver::GetCurrentVelocity()
{
	return m.flowsolver->GetCurrentVelocity();
}

mfem::Mesh* MFEMSolver::GetMFEMMesh()
{
	return m.mesh;
}

mfem::Mesh* MFEMSolver::GetMFEMfluidmesh()
{
	return m.fluid_mesh;
}

mfem::Mesh* MFEMSolver::GetMFEMsolidmesh()
{
	return m.solid_mesh;
}
#endif

bool MFEMSolver::ThreeFieldFormulation()
{
   return m.three_field;
}

bool MFEMSolver::ImmersedSolid()
{
   return m.immersed_solid;
}

int MFEMSolver::SolidSolverType()
{
   return m.solidsolvertype;
}

void MFEMSolver::AddVelDirichletBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr)
{
	m.flowsolver->AddVelDirichletBC(coeff, attr);
}

void MFEMSolver::AddVelDirichletBC(mfem::Coefficient* coeff, mfem::Array<int>& attr, int component)
{
	m.flowsolver->AddVelDirichletBC(coeff, attr, component);
}

void MFEMSolver::AddNoSlipBC(mfem::Array<int>& attr)
{
   m.flowsolver->AddNoSlipBC(attr);
}

void MFEMSolver::SetNodalPressure(FEParamDouble pval, mfem::Vector &node_pos)
{
   if (pval.isConst() == false)
   {
      feLogError("Pressure value must be constant");
      assert(false);
   }
   double pressure_value = pval.constValue();
   m.flowsolver->AddPressureNodeCondition(pressure_value, node_pos);
}

void MFEMSolver::AddPrescribedVelocityComponentBC(mfem::Coefficient* coeff, mfem::Array<int>& attr, int dof)
{
   m.flowsolver->AddVelDirichletBC(coeff, attr, dof);
}

void MFEMSolver::AddTractionBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr)
{
	m.flowsolver->AddTractionBC(coeff, attr);
}

void MFEMSolver::AddTractionFreeBC(mfem::Array<int>& attr)
{
	m.flowsolver->AddTractionFreeBC(attr);
} 

void MFEMSolver::AddPressureBC(mfem::Coefficient* coeff, mfem::Array<int>& attr)
{
   m.flowsolver->AddPressureBC(coeff, attr);
}

void MFEMSolver::AddPeriodicBC(mfem::Array<int>& attr)
{
   m.flowsolver->AddPeriodicBC(attr);
}

void MFEMSolver::AddBackflowStabilizationBC(FEParamDouble beta, mfem::Array<int>& attr)
{
   if (beta.isConst() == false)
   {
      feLogError("Backflow stabilization coefficient must be constant");
      assert(false);
   }
   double beta_val = beta.constValue();
   m.flowsolver->UseBackflowStabilization(attr,beta_val);
}

void MFEMSolver::AddResistanceBC(FEParamDouble R, FEParamDouble p0, mfem::Array<int>& attr)
{
   if (R.isConst() == false)
   {
      feLogError("Resistance coefficient must be constant");
      assert(false);
   }
   if (p0.isConst() == false)
   {
      feLogError("Initial pressure must be constant");
      assert(false);
   }
   double R_val = R.constValue();
   double p0_val = p0.constValue();
   m.flowsolver->AddResistanceBC(attr,R_val,p0_val);
}

void MFEMSolver::AddWindkesselBC(FEParamDouble R, FEParamDouble Rd, FEParamDouble C, FEParamDouble p0, FEParamDouble pd, mfem::Array<int>& attr)
{
   if (R.isConst() == false)
   {
      feLogError("Proximal resistance coefficient must be constant");
      assert(false);
   }
   if (Rd.isConst() == false)
   {
      feLogError("Distal resistance coefficient must be constant");
      assert(false);
   }
   if (C.isConst() == false)
   {
      feLogError("Compliance coefficient must be constant");
      assert(false);
   }
   if (p0.isConst() == false)
   {
      feLogError("Initial pressure must be constant");
      assert(false);
   }
   if (pd.isConst() == false)
   {
      feLogError("Distal pressure must be constant");
      assert(false);
   }
   double R_val = R.constValue();
   double Rd_val = Rd.constValue();
   double C_val = C.constValue();
   double p0_val = p0.constValue();
   double pd_val = pd.constValue();

   // compute initial flow rate
   #ifdef MFEM_USE_MPI
   ParGridFunction *u_ic = m.flowsolver->GetCurrentVelocity();
   #else
   GridFunction *u_ic = m.flowsolver->GetCurrentVelocity();
   #endif
   double Q0 = m.flowsolver->ComputeMassFlowThruSurf(*u_ic, attr, false, m.fluid_density, false);

   m.flowsolver->AddWindkesselBC(attr,R_val,Rd_val,C_val,p0_val,pd_val,Q0);
}

void MFEMSolver::AddBodyForce(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr)
{
   m.flowsolver->AddBodyForce(coeff, attr);
}

double MFEMSolver::FluidDensity() const
{
	return m.fluid_density;
}

double MFEMSolver::KinematicViscosity() const
{
	return m.kinvis;
}

double MFEMSolver::DynamicViscosity() const
{
	return m.dyn_viscosity;
}

bool MFEMSolver::Init()
{
	// create the mfem mesh from the FEBio mesh
	FEModel *fem = GetFEModel();

   // Redirect mfem::out to the custom stream buffer that writes to the log file
   // TODO: Cleanup the output to the log file
   m.customBuffer.originalBuffer = mfem::out.rdbuf();
   m.customBuffer.fe_model = fem;
   std::ostream customStream(&m.customBuffer);
   mfem::out.rdbuf(customStream.rdbuf());

   // create the mesh
	FEMesh& femesh = fem->GetMesh();
	m.mesh = CreateMFEMMesh3D(femesh);
   // m.mesh->EnsureNodes();
	if (m.mesh == nullptr) return false;
   if (m.immersed_solid)
   {
      // create fluid and solid submeshes from mfem version of FEBio mesh
      bool mfem_solid_material = false;
      bool febio_solid_material = false;
      Array<int> solid_domain_attr(m.mesh->attributes.Max());
      solid_domain_attr = 0;
      for (int i = 0; i < femesh.Domains(); ++i)
      {
         FEDomain& dom = femesh.Domain(i);
         int dom_attr = m.mesh->GetDomainAttribute(&dom);
         if (dom.IsActive() && dynamic_cast<MFEMSolidMaterial*>(dom.GetMaterial()))
         {
            // solid_domain_attr[i] = i+1;
            solid_domain_attr[i] = dom_attr;
            mfem_solid_material = true;
         }
         if (dom.IsActive() && dynamic_cast<FEElasticMaterial*>(dom.GetMaterial()))
         {
            // solid_domain_attr[i] = i+1;
            solid_domain_attr[i] = dom_attr;
            febio_solid_material = true;
         }
         if(dom.IsActive() && dynamic_cast<MFEMNewtonianFluid*>(dom.GetMaterial()))
         {
            Array<int> domain_attr(m.mesh->attributes.Max());
            domain_attr = 0;
            // domain_attr[i] = i+1;
            domain_attr[i] = dom_attr;
            m.fluid_mesh = new SubMesh(SubMesh::CreateFromDomain(*m.mesh,domain_attr));
            if (m.fluid_mesh == nullptr) return false;
         }
      }

      // create solid mesh as sum of all solid domains, this will preserve solid domain numbering
      if (mfem_solid_material && febio_solid_material)
      {
         feLogError("Cannot have both FEBio and MFEM solid materials in the same mesh");
         return false;
      }
      else
      {
         m.solid_mesh = new SubMesh(SubMesh::CreateFromDomain(*m.mesh,solid_domain_attr));
         if (m.solid_mesh == nullptr) return false;
      }
   }

   // print mesh details
   if (m.immersed_solid)
   {
      feLog("-------------------------------\nFluid Mesh\n");
      m.fluid_mesh->PrintCharacteristics();
      feLog("-------------------------------\nSolid Mesh\n");
      m.solid_mesh->PrintCharacteristics();
      feLog("-------------------------------\n");
   }
   else
   {
      feLog("-------------------------------\nFluid Mesh\n");
      m.mesh->PrintCharacteristics();
      feLog("-------------------------------\n");
   }
   
   #ifdef MFEM_USE_MPI
   // Partition the meshes
   auto fluid_serial_mesh = m.immersed_solid ? static_cast<Mesh*>(m.fluid_mesh) : static_cast<Mesh*>(m.mesh);
   m.fluid_pmesh = new ParMesh(MPI_COMM_WORLD, *fluid_serial_mesh);
   if (m.immersed_solid) m.solid_pmesh = new ParMesh(MPI_COMM_WORLD, *m.solid_mesh);
   
   if (m.immersed_solid && m.solidsolvertype == SolidSolverType::FEBio)
   {
      // create ParFEModel
      m.par_fem = new ParFEModel(fem);
      m.par_fem->SetMFEMSerialMesh(m.mesh);
      m.par_fem->SetMFEMSubMesh(m.solid_mesh);
      m.par_fem->SetupBCs(); // Probably not needed if only clear mesh and not surfaces associated with each bc!
   }
   #endif

   // get fluid material
   for (int i=0; i<fem->Materials(); ++i)
   {
      FEMaterial* mat = fem->GetMaterial(i);
      if (dynamic_cast<MFEMNewtonianFluid*>(mat))
      {
         m.fluid_density = dynamic_cast<MFEMNewtonianFluid*>(mat)->Density();
         m.dyn_viscosity = dynamic_cast<MFEMNewtonianFluid*>(mat)->DynamicViscosity();
         m.kinvis = dynamic_cast<MFEMNewtonianFluid*>(mat)->KinematicViscosity();
         break;
      }
   } 

	// create the mfem solver class
   if (!m.stabilized_VMS)
   {
      if (m.orderv <= 1.0) 
      {
         feLogError("Mixed order method requires velocity order greater than or equal to 2");
         return false;
      }
      if (m.orderp != m.orderv - 1.0)
      {
         feLogError("Mixed order method requires pressure order to be one less than velocity order");
         return false;
      }
   }
   
   if (m.linsolvertype == m.LinearSolverType::PARDISO)
   {
      #ifdef MFEM_USE_MPI
      feLogError("PARDISO solver is not supported in parallel");
      return false;
      #else
      m.flowsolver = new FEBio_INS_Solver(fem, m.immersed_solid ? static_cast<Mesh*>(m.fluid_mesh) : static_cast<Mesh*>(m.mesh), (int) m.orderv, (int) m.orderp, m.kinvis, m.fluid_density, m.orderingtype, m.immersed_solid);
      #endif
   }
   else
   {
      #ifdef MFEM_USE_MPI
      m.flowsolver = new INS_Solver(m.fluid_pmesh, (int) m.orderv, (int) m.orderp, m.kinvis, m.fluid_density, m.orderingtype, m.immersed_solid);
      #else
      m.flowsolver = new INS_Solver(m.immersed_solid ? static_cast<Mesh*>(m.fluid_mesh) : static_cast<Mesh*>(m.mesh), (int) m.orderv, (int) m.orderp, m.kinvis, m.fluid_density, m.orderingtype, m.immersed_solid);
      #endif
   }
   
	// set rho_inf for time stepping method
	m.flowsolver->SetRhoInf(m.rhoinf);

	// set tolerances
	m.flowsolver->SetLinearSolveAbsTol(m.lin_atol);
	m.flowsolver->SetLinearSolveRelTol(m.lin_rtol);
	m.flowsolver->SetLinearSolveMaxIt(m.lin_maxIt);

	m.flowsolver->SetNewtonSolveAbsTol(m.atol);
	m.flowsolver->SetNewtonSolveRelTol(m.rtol);
	m.flowsolver->SetNewtonSolveMinIt(m.minIt);
	m.flowsolver->SetNewtonSolveMaxIt(m.maxIt);

	// set linear solver and preconditioner
   if (m.linsolvertype == m.LinearSolverType::PARDISO || m.linsolvertype == m.LinearSolverType::UMFPACK || m.linsolvertype == m.LinearSolverType::CPARDISO)
   {
      #ifdef MFEM_USE_MPI
      MFEM_ASSERT(m.linsolvertype == m.LinearSolverType::CPARDISO, "CPARDISO is the only supported direct linear solver in parallel");
      #else
      MFEM_ASSERT(m.linsolvertype != m.LinearSolverType::CPARDISO, "CPARDISO is not supported in serial");
      #endif

      m.direct_solver = true;
      m.flowsolver->UseDirectSolver();
   }
   else 
   {
      m.direct_solver_prec = m.prec_inner_solver_type == m.PrecInnerSolverType::Direct ? true : false;
      m.flowsolver->SetPreconditioner(m.prec_type, m.direct_solver_prec);
   }

	// set nonlinear form type
	m.flowsolver->SetNonlinearFormType(m.nonlinear_form);

   // use total integrator
	if (m.total_integ) m.flowsolver->UseTotalIntegrator();

	// stabilizations
   // TODO: setup options for each, for now just VMS
   if (m.stabilized_VMS)
   {
      switch (m.stab_type) {
       case ins_solver::StabilizationType::NONE:
            break;
		 case ins_solver::StabilizationType::SUPG:
            if (m.orderv != m.orderp - 1) 
            {
                  feLogError("using SUPG only requires mixed order elements");
                  return false;
            }
            m.flowsolver->UseSUPG();
            break;

		 case ins_solver::StabilizationType::SUPG_PSPG:
            if (m.orderv != m.orderp) 
            {
                  feLogError("SUPG/PSPG requires equal order elements");
                  return false;
            }
            m.flowsolver->UseSUPG();
            m.flowsolver->UsePSPG();
            break;

		 case ins_solver::StabilizationType::VMS:
         default:
            if (m.orderv != m.orderp) 
            {
                  feLogError("VMS requires equal order elements");
                  return false;
            }
            m.flowsolver->UseVMS();
            break;
      }
   }
	if (m.stabilized_LSIC) 
   {
      // turned on by default with VMS
      m.flowsolver->UseMassStabilization();
   }

	// set viscous term type, default is StressForm
	if (m.LaplaceForm)
	{
		m.flowsolver->SetViscousTermType(ins_solver::ViscousTermType::LaplaceForm);
	}

   // TODO: In parallel these need to be set after creating the ParFESpaces
   // so the correct boundary dofs get selected on each processor!

   // Here we assume the immersed solid has no initial or boundary conditions set
   // Set the initial conditions.
	for (int i = 0; i < fem->InitialConditions(); ++i)
	{
		MFEMInitialCondition* ic = dynamic_cast<MFEMInitialCondition*>(fem->InitialCondition(i));
		if (ic && ic->IsActive())
		{
			ic->Apply(this);
		}
	}

	// set boundary conditions and loads
	for (int i=0; i<fem->BoundaryConditions(); ++i)
	{
		MFEMBoundaryCondition* bc = dynamic_cast<MFEMBoundaryCondition*>(fem->BoundaryCondition(i));
		if (bc && bc->IsActive())
		{
			FESurface* surf = bc->GetSurface();
         int id = 0;
         if (surf)
         {
            FEFacetSet* set = surf->GetFacetSet();
            id = m.mesh->GetSurfaceAttribute(set); assert(id > 0);
         }
         else
         {
            id = bc->GetBoundaryID(); assert(id > 0);
         }		
         
			bc->Apply(this, id - 1);
		}
	}
	for (int i = 0; i < fem->ModelLoads(); ++i)
	{
      // nodal loads
      MFEMNodalLoad* nl = dynamic_cast<MFEMNodalLoad*>(fem->ModelLoad(i));
      if (nl && nl->IsActive())
      {
         nl->Apply(this);
      }

      // surface loads
		MFEMSurfaceLoad* sl = dynamic_cast<MFEMSurfaceLoad*>(fem->ModelLoad(i));
		if (sl && sl->IsActive())
		{
         FESurface* surf = &sl->GetSurface();
         int id = 0;
         if (surf)
         {
            FEFacetSet* set = surf->GetFacetSet();
            id = m.mesh->GetSurfaceAttribute(set); assert(id > 0);
         }
         else
         {
            id = sl->GetBoundaryID(); assert(id > 0);
         }

         sl->Apply(this, id - 1);
		}

      // body loads
      MFEMBodyLoad* bl = dynamic_cast<MFEMBodyLoad*>(fem->ModelLoad(i));
      if (bl && bl->IsActive())
      {
         std::vector<int> &dom_id_list = bl->GetDomainIDList();
         if (dom_id_list.size() != 0)
         {
            for (int j=0; j<dom_id_list.size(); ++j)
            {
               int id = dom_id_list[j]; assert(id > 0);
               bl->Apply(this, id - 1);
            }
         }
         else
         {
            for (int j=0; j<bl->Domains(); ++j)
            {
               FEDomain* dom = bl->Domain(j);
               int id = m.mesh->GetDomainAttribute(dom); assert(id > 0);
               bl->Apply(this, id - 1);
            }
         }                  
      }
	}

   // Zero mean pressure
   if (m.zero_mean_pressure) m.flowsolver->AddZeroMeanPressureCondition();

   // Fluid grid functions
	m.u_gf = m.flowsolver->GetCurrentVelocity();
	m.p_gf = m.flowsolver->GetCurrentPressure();

   // Add immersed solid
   if (m.immersed_solid)
   {
      MFEM_ASSERT(fem->Materials() < 3, "Immersed solid solver only supports one material at this time");

      // create solid material
      MFEMSolidMaterial* solid_mat = nullptr;
      FEElasticMaterial* febio_elast_mat = nullptr;
      for (int i=0; i<fem->Materials(); ++i)
      {
         FEMaterial* fe_mat = fem->GetMaterial(i);
         
         // MFEM material
         solid_mat = dynamic_cast<MFEMSolidMaterial*>(fe_mat);
         if (solid_mat) { break; }
         
         // FEBio material
         febio_elast_mat = dynamic_cast<FEElasticMaterial*>(fe_mat);
         if (febio_elast_mat)
         {
            MFEM_ASSERT(m.form_type == ImmersedSolid::form_type::UpdatedLagrangian, "FEBio materials only work with Updated Lagrangian");
            break;          
         }
      }
      if (!solid_mat && !febio_elast_mat)
      {
         feLogError("No or invalid immersed solid material found");
         return false;
      }

      // determine which material
      FEScalarCoefficient *init_rho_coeff = nullptr;
      if (solid_mat)
      {
         // MFEM Material
         init_rho_coeff = solid_mat->GetDensityCoeff();
         if (typeid(*solid_mat) == typeid(MFEMCompressibleNeoHookeanMaterial))
         {
            MFEMCompressibleNeoHookeanMaterial* cnh_mat = static_cast<MFEMCompressibleNeoHookeanMaterial*>(solid_mat);
            m.mat_type = new CompressibleNeoHookean(m.form_type,m.solid_mesh->Dimension(),*cnh_mat->GetC1Coeff(),*cnh_mat->GetC0Coeff());
         }
         else if (typeid(*solid_mat) == typeid(MFEMNearlyIncompressibleNeoHookeanMaterial))
         {
            MFEMNearlyIncompressibleNeoHookeanMaterial* nih_mat = static_cast<MFEMNearlyIncompressibleNeoHookeanMaterial*>(solid_mat);
            m.mat_type = new NearlyIncompressibleNeoHookean(m.form_type,m.solid_mesh->Dimension(),m.vol_type,*nih_mat->GetBulkModulus(),*nih_mat->GetC0Coeff());
         }
         else if (typeid(*solid_mat) == typeid(MFEMNearlyIncompressibleAnnularDisk))
         {
            MFEMNearlyIncompressibleAnnularDisk* niad_mat = static_cast<MFEMNearlyIncompressibleAnnularDisk*>(solid_mat);
            double disk_w = niad_mat->GetRadialThickness(); // constant thickness
            m.mat_type = new NearlyIncompressibleAnnularDisk(m.form_type,m.solid_mesh->Dimension(),m.vol_type,*niad_mat->GetBulkModulus(),*niad_mat->GetC0Coeff(),disk_w);
         }
         else
         {
            feLogError("Unsupported solid material type");
            return false;
         }
      }
      else if (febio_elast_mat)
      {
         // FEBio material
         m.febio_material = true;
         
         if (m.solidsolvertype == SolidSolverType::MFEM)
         {
            // MFEM Solver with FEBio elastic material
            FEBioElasticMaterial* feb_mat = new FEBioElasticMaterial(febio_elast_mat);
            init_rho_coeff = feb_mat->GetDensityCoeff();
            m.mat_type = feb_mat;
         }
      }
      else
      {
         feLogError("No immersed solid material found");
         return false;
      }      

      // add immersed solid to the flow solver
      if (m.solidsolvertype == SolidSolverType::MFEM)
      {
         MFEM_ASSERT(fem->Materials() < 3, "MFEM solid solver only supports one material");

         #ifdef MFEM_USE_MPI
         m.flowsolver->AddImmersedSolid(m.solid_pmesh,m.orderd,init_rho_coeff,m.mat_type,m.form_type,m.subtract_viscous_stress,m.stabs_scale);
         #else
         m.flowsolver->AddImmersedSolid(m.solid_mesh,m.orderd,init_rho_coeff,m.mat_type,m.form_type,m.subtract_viscous_stress,m.stabs_scale);
         #endif

         // immersed solid grid functions
         m.d_gf = m.flowsolver->GetImmersedSolidDisplacement();
         m.vib_gf = m.flowsolver->GetImmersedSolidVelocity();
      }
      else if (m.solidsolvertype == SolidSolverType::FEBio)
      {
         MFEM_ASSERT(m.febio_material, "Immersed FEBio solid solver only supports FEBio materials");

         auto body_force_terms = m.flowsolver->GetBodyForces();
         #ifdef MFEM_USE_MPI
         auto fluid_vfes = m.flowsolver->GetCurrentVelocity()->ParFESpace();
         m.immersed_febio_solid = new ImmersedFEBioSolid(fem,m.solid_pmesh,m.orderd,init_rho_coeff,m.fluid_density,m.dyn_viscosity,body_force_terms,fluid_vfes,m.fluid_pmesh,m.mat_type,m.orderingtype,m.form_type,m.subtract_viscous_stress);
         #else
         auto fluid_vfes = m.flowsolver->GetCurrentVelocity()->FESpace();
         m.immersed_febio_solid = new ImmersedFEBioSolid(fem,m.solid_mesh,m.orderd,init_rho_coeff,m.fluid_density,m.dyn_viscosity,body_force_terms,fluid_vfes,m.fluid_mesh,m.mat_type,m.orderingtype,m.form_type,m.subtract_viscous_stress);
         #endif

         // add immersed solid to the flow solver
         m.flowsolver->AddImmersedSolid(m.immersed_febio_solid,m.stabs_scale);

         // immersed solid grid functions
         m.d_gf = m.flowsolver->GetImmersedSolidDisplacement();
         m.vib_gf = m.flowsolver->GetImmersedSolidVelocity();
      }
      else
      {
         feLogError("Invalid immersed solid solver type");
         return false;
      }         
   }   

   // create transfer maps and global functions
   if(m.immersed_solid)
   {
      #ifndef MFEM_USE_MPI
      int vdim = m.flowsolver->GetCurrentVelocity()->FESpace()->GetVDim();
      int ordering = m.flowsolver->GetCurrentVelocity()->FESpace()->GetOrdering();
      auto vfec = m.flowsolver->GetCurrentVelocity()->FESpace()->FEColl();
      auto pfec = m.flowsolver->GetCurrentPressure()->FESpace()->FEColl();
      m.vfespace_global = new FiniteElementSpace(static_cast<Mesh*>(m.mesh),vfec,vdim,ordering);
      m.pfespace_global = new FiniteElementSpace(static_cast<Mesh*>(m.mesh),pfec,1,ordering);
   
      m.u_gf_global.SetSpace(m.vfespace_global);
      m.u_gf_global = 0.0;
      m.p_gf_global.SetSpace(m.pfespace_global);
      m.p_gf_global = 0.0;
      m.d_gf_global.SetSpace(m.vfespace_global);
      m.d_gf_global = 0.0;
      m.vib_gf_global.SetSpace(m.vfespace_global);
      m.vib_gf_global = 0.0;

      m.submeshfluid_to_global_v = new TransferMap(*m.u_gf,m.u_gf_global);
      m.submeshfluid_to_global_p = new TransferMap(*m.p_gf,m.p_gf_global);
      m.submeshsolid_to_global_d = new TransferMap(*m.d_gf,m.d_gf_global);
      m.global_to_submeshsolid_d = new TransferMap(m.d_gf_global,*m.d_gf);
      
      if (m.solidsolvertype == SolidSolverType::FEBio)
      {
         // set submesh mapping data
         m.immersed_febio_solid->SetSubMeshMaps(m.submeshsolid_to_global_d,m.global_to_submeshsolid_d,&m.vib_gf_global);
      }
      #endif
   }

   #ifdef MFEM_USE_MPI
   // build local FEBio solid Mesh from mfem solid ParMesh
   femesh.Clear(); // clear FEBio serial mesh
   if (m.immersed_solid && m.solidsolvertype == SolidSolverType::FEBio)
   {
      ParFiniteElementSpace *dfes = m.flowsolver->GetImmersedSolidDisplacement()->ParFESpace();
      m.par_fem->BuildLocalFEBioMesh(dfes);
   }
   
   // delete the serial MFEM meshes
   delete m.mesh;
   delete m.fluid_mesh;
   delete m.solid_mesh;  
   #endif

	// setup output in paraview format
   std::string fname = "results_mfem_solver";
   #ifdef MFEM_USE_MPI
   fname += "_np" + std::to_string(mfem::Mpi::WorldSize());
   #endif
   if (m.immersed_solid)
   {
      std::string fluid_filename = fname + "/fluid";
      #ifdef MFEM_USE_MPI
      m.pvdc = new ParaViewDataCollection(fluid_filename, m.fluid_pmesh);
      #else
      m.pvdc = new ParaViewDataCollection(fluid_filename, m.fluid_mesh);
      #endif      
   }
   else
   {
      #ifdef MFEM_USE_MPI
      m.pvdc = new ParaViewDataCollection(fname, m.fluid_pmesh);
      #else
      m.pvdc = new ParaViewDataCollection(fname, m.mesh);
      #endif      
   }	
	m.pvdc->SetDataFormat(VTKFormat::BINARY32);
	m.pvdc->SetHighOrderOutput(true);
	m.pvdc->SetLevelsOfDetail(m.orderv);
	m.pvdc->SetCycle(0);
	m.pvdc->SetTime(0);
	m.pvdc->RegisterField("Velocity", m.u_gf);
	m.pvdc->RegisterField("Pressure", m.p_gf);
   if(m.immersed_solid)
   {
      if (m.flowsolver->GetImmersedSolidIndicator()->FESpace() != nullptr)
      {
         #ifdef MFEM_USE_MPI
         ParGridFunction *s_gf = m.flowsolver->GetImmersedSolidIndicator();
         #else
         GridFunction *s_gf = m.flowsolver->GetImmersedSolidIndicator();
         #endif
         m.pvdc->RegisterField("Indicator", s_gf);
      }      
   }
   // TODO: add more outputs (curlv, wss, others)
	m.pvdc->Save();

   if (m.immersed_solid)
   {
      std::string solid_filename = fname + "/solid";
      #ifdef MFEM_USE_MPI
      m.sdc = new ParaViewDataCollection(solid_filename, m.solid_pmesh);
      #else
      m.sdc = new ParaViewDataCollection(solid_filename, m.solid_mesh);  
      #endif
      m.sdc->SetDataFormat(VTKFormat::BINARY);
      m.sdc->SetHighOrderOutput(true);
      m.sdc->SetLevelsOfDetail(m.orderd);
      m.sdc->RegisterField("Displacement", m.d_gf);
      m.sdc->RegisterField("Velocity", m.vib_gf);
      m.sdc->SetCycle(0);
      m.sdc->SetTime(0);
      m.sdc->Save();
   }

	double dt = fem->GetTime().timeIncrement;
	m.flowsolver->Setup(dt);
   
	return true;
}

bool MFEMSolver::SolveStep()
{
	// Get the FEBio model and mesh objects
	FEModel& fem = *GetFEModel();
	FEMesh& femesh = fem.GetMesh();

	FEAnalysis* analysis = fem.GetCurrentStep();
	int timeStep = analysis->m_ntimesteps + 1; // current time step number

	double t = fem.GetTime().currentTime; // t = t_n+1
	double dt = fem.GetTime().timeIncrement;

   // solve a time step
	m.flowsolver->Step(t, dt, timeStep);
   
	// save data
   if ((timeStep) % m.step_save == 0)
   {
      m.pvdc->SetCycle(timeStep);
      m.pvdc->SetTime(t);
      m.pvdc->Save();
      if (m.immersed_solid)
      {
         m.sdc->SetCycle(timeStep);
         m.sdc->SetTime(t);
         m.sdc->Save();
      }
   }	

   // TODO how to store in parallel?
   #ifndef MFEM_USE_MPI

	// copy data to FEBio mesh
   if(m.immersed_solid)
   {
      // transfer grid functions from submesh to parent febio mesh
      m.u_gf_global = 0.0;
      m.p_gf_global = 0.0;
      m.d_gf_global = 0.0;
      m.vib_gf_global = 0.0;
      m.submeshfluid_to_global_v->Transfer(*m.u_gf, m.u_gf_global);
      m.submeshfluid_to_global_p->Transfer(*m.p_gf, m.p_gf_global);
      m.submeshsolid_to_global_d->Transfer(*m.d_gf, m.d_gf_global);
      m.submeshsolid_to_global_d->Transfer(*m.vib_gf, m.vib_gf_global);

      Vector& p = m.p_gf_global;
      Vector& v = m.u_gf_global;
      int dofP = fem.GetDOFIndex("p");
      int dofVX = fem.GetDOFIndex("vx");
      int dofVY = fem.GetDOFIndex("vy");
      int dofVZ = fem.GetDOFIndex("vz");
      int dofUX, dofUY, dofUZ, dofVIMX, dofVIMY, dofVIMZ;
      Vector& d = m.d_gf_global;
      Vector& vib = m.vib_gf_global;
      dofUX = fem.GetDOFIndex("im_ux");
      dofUY = fem.GetDOFIndex("im_uy");
      dofUZ = fem.GetDOFIndex("im_uz");
      dofVIMX = fem.GetDOFIndex("im_vx");
      dofVIMY = fem.GetDOFIndex("im_vy");
      dofVIMZ = fem.GetDOFIndex("im_vz");

      int NN = femesh.Nodes();
      for (int i = 0; i < NN; ++i)
      {
         FENode& node = femesh.Node(i);
         node.set(dofP, p[i]);

         // for vector unknowns, check dof ordering type
         int idx = m.orderingtype == Ordering::byNODES ? i + 0*NN : i*3 + 0;
         int idy = m.orderingtype == Ordering::byNODES ? i + 1*NN : i*3 + 1;
         int idz = m.orderingtype == Ordering::byNODES ? i + 2*NN : i*3 + 2;
         node.set(dofVX, v[idx]);
         node.set(dofVY, v[idy]);
         node.set(dofVZ, v[idz]);

         node.set(dofUX, d[idx]);
         node.set(dofUY, d[idy]);
         node.set(dofUZ, d[idz]);
         node.set(dofVIMX, vib[idx]);
         node.set(dofVIMY, vib[idy]);
         node.set(dofVIMZ, vib[idz]);
      }
      fem.Update();
   }
   else
   {
      Vector& p = *m.p_gf;
      Vector& v = *m.u_gf;
      int dofP = fem.GetDOFIndex("p");
      int dofVX = fem.GetDOFIndex("vx");
      int dofVY = fem.GetDOFIndex("vy");
      int dofVZ = fem.GetDOFIndex("vz");
      int NN = femesh.Nodes();
      for (int i = 0; i < NN; ++i)
      {
         FENode& node = femesh.Node(i);
         node.set(dofP, p[i]);
         
         // for vector unknowns, check dof ordering type
         int idx = m.orderingtype == Ordering::byNODES ? i + 0*NN : i*3 + 0;
         int idy = m.orderingtype == Ordering::byNODES ? i + 1*NN : i*3 + 1;
         int idz = m.orderingtype == Ordering::byNODES ? i + 2*NN : i*3 + 2;
         node.set(dofVX, v[idx]);
         node.set(dofVY, v[idy]);
         node.set(dofVZ, v[idz]);
      }
      fem.Update();
   }
   #endif

	// print out mfem timing data
	if (timeStep == analysis->m_ntime)
	{
		m.flowsolver->PrintTimingData();
	}

	return true;
}
