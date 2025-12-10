#include "MFEMSolver.h"

#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FEBioLib/FEBioModel.h>

#include "MFEMMesh.h"
#include "MFEMInitialCondition.h"
#include "MFEMBoundaryCondition.h"
#include "MFEMLoad.h"
#include "MFEMMaterial.h"
#include "../MFEMLinAlg/FEBioLinearSolver.h"
#include "../MFEMImSolid/MFEMImmersedElasticDomain.h"
#include "../MFEMImSolid/MFEMImmersedFEBioSolid.h"
#include "../MFEMUtilities/MFEMOutputTimeSeries.h"
#include "../MFEMUtilities/MFEMParFEMeshBuilder.h"

#include "../MFEMFSI/Navier-Stokes/IncompressibleFlowSolver.hpp"
#include "../MFEMFSI/Navier-Stokes/ifs_solver_enum.hpp"
#include "../MFEMFSI/ImmersedFSI/ImmersedSolid.hpp"
#include "../MFEMFSI/ImmersedFSI/Materials/materials.hpp"

#include <typeinfo>
#include <filesystem>
#include <sstream>
#include <iomanip>

using namespace ifs_solver;
using namespace std;
using namespace std::placeholders;
using namespace immersed_solid;
using namespace mat;
using namespace mfem;

// additional functions
double ComputeVectorH1Error(GridFunction *u_gf, VectorCoefficient *ex_u, MatrixCoefficient *ex_gradu, const IntegrationRule *irs[]);   
void ComputeQCriterion(GridFunction &u_gf, GridFunction &q_gf);    
void ComputeAverageVelocity(int form_type, GridFunction &u_gf, GridFunction &v_gf, Vector &avg_velo, double &volume);
template<typename T>
T sq(T x)
{
   return x * x;
};  

#ifndef MFEM_USE_MPI
// we override this so we can use the FEBio linear solvers
class FEBio_INS_Solver : public IncompressibleFlowSolver
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
		bool immersed_flag = false) : IncompressibleFlowSolver(mesh, orderv, orderp, kin_vis, fluid_density_, ordering_type_, immersed_flag) 
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
	IncompressibleFlowSolver* flowsolver = nullptr; // MFEM flow solver
   #ifdef MFEM_USE_MPI
   ParGridFunction* p_gf = nullptr;
   ParGridFunction* u_gf = nullptr;
   #else
	GridFunction* p_gf = nullptr; // velocity grid function
	GridFunction* u_gf = nullptr; // pressure grid function
   #endif
	ParaViewDataCollectionBase* pvdc = nullptr; // paraview data collection
   #ifdef MFEM_USE_MPI
   ParMesh* fluid_pmesh = nullptr;
   #endif
   //additional output options
   bool wss = false;
   bool vorticity = false;
   bool q_criterion = false;
   bool output_surface_pressures = false;
   std::unique_ptr<SurfacePressureTimeSeries> surf_press_ts;
   bool output_surface_flowrates = false;
   std::unique_ptr<SurfaceFlowRateTimeSeries> surf_flow_ts;
   struct SurfaceData
   {
      std::string name;
      int attribute;
   };
   mutable std::vector<SurfaceData> output_surfaces;
   #ifdef MFEM_USE_MPI
   ParGridFunction* wss_gf = nullptr; // wall shear stress grid function
   ParGridFunction vorticity_gf; // vorticity grid function
   ParGridFunction q_criterion_gf; // Q-criterion grid function
   #else
   GridFunction* wss_gf = nullptr; // wall shear stress grid function
   GridFunction vorticity_gf; // vorticity grid function
   GridFunction q_criterion_gf; // Q-criterion grid function
   #endif

   // immersed solid
   SubMesh* solid_mesh = nullptr; // immersed solid mesh
   #ifdef MFEM_USE_MPI
   ParGridFunction *d_gf = nullptr; // immersed solid displacement grid function
   ParGridFunction *vib_gf = nullptr; // immersed solid velocity grid function
   ParGridFunction *a_gf = nullptr; // immersed solid acceleration grid function
   #else
   GridFunction *d_gf = nullptr; // immersed solid displacement grid function
   GridFunction *vib_gf = nullptr; // immersed solid velocity grid function
   GridFunction *a_gf = nullptr; // immersed solid acceleration grid function
   #endif
   ParaViewDataCollectionBase *sdc = nullptr; // immersed solid paraview data collection
   int solidsolvertype = 1; // 0 = MFEM, 1 = FEBio (default is MFEM)
   ImmersedFEBioSolid *immersed_febio_solid = nullptr; // not owned
   bool three_field = false; // flag for 3-field formulation
   #ifdef MFEM_USE_MPI
   ParMesh* solid_pmesh = nullptr;
   #endif
   std::map<FESolidDomain*, int> dom_attr_map; // map from immersed solid domain to mfem attribute
   GridFunction *vonMisesStress_gf = nullptr; // von Mises stress grid function

   // Output
   std::string fname = "results_mfem_solver";
   mfem::VTKFormat output_format = mfem::VTKFormat::BINARY32; // ASCII, BINARY, or BINARY32
   int compression_level = -1; // zlib compression -1 default, max 9

   // Restarting
   std::string restart_folder = "RestartFiles";
   bool write_restarts = false;
   int restart_step_save = 1;
   bool restarting = false;
   int restart_step = 1;
   double restart_t = 0.0;
   double restart_dt = 0.0;

   // Starting at a specific time and step
   bool start_from_specific_time = false;
   double start_sim_time = 0.0;
   int start_sim_step = 0;

   // Write final state as IC for next run
   bool write_final_state_as_IC = false;
   std::string final_state_folder = "InitialConditions";

   // Pre-partitioning of the mesh
   // Large meshes can be pre-partitioned, then only the local piece is read on each rank
   // rather than the entire mesh read on each rank then performing the partitioning
   bool pre_partition = false; // flag for pre-partitioning the mesh
   int partition_opt = 0; // 0 = PartitionMesh, 1 = ReadParallelMesh
   int np = 1; // number of ranks for pre-partitioning
   enum ParallelMesh { PartitionMesh = 0, ReadParallelMesh = 1 };

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
   int stab_type = ifs_solver::StabilizationType::VMS;
   int orderingtype= Ordering::byVDIM;

   int nonlinear_form = ifs_solver::NonlinearFormType::Convective;
   int viscous_term_type = ifs_solver::ViscousTermType::StressForm;
   bool LaplaceForm = false;
   bool total_integ = true;
   bool zero_mean_pressure = false;

   // Time Stepping Options
   int ts_method = ifs_solver::TimeSteppingMethod::GeneralizedAlpha;
   double rhoinf = 0.5;
   int step_save = 1; // frequency of storing data
   int steady_state_step_save = 1; // frequency of storing data after steady state
   double steady_state_time = 0.0;
   bool steady_state_save = false;

   // Nonlinear Solver Options
   int minIt = 3;
   int maxIt = 10;
   double atol = 1.0e-10;
   double rtol = 1.0e-3;
   bool ns_no_error_nonconvergence = false;

   // Linear Solver Options
   // int lin_minIt;
   int lin_maxIt = 200;
   double lin_atol = 1.0e-10;
   double lin_rtol = 1.0e-3;
   bool direct_solver = false;
   // PARDISO = OpenMP (serial code only), CPARDISO = MPI (parallel code only)
   enum LinearSolverType { UMFPACK = 0, PARDISO = 1, FGMRES = 2, CPARDISO = 3 }; 
   int linsolvertype = LinearSolverType::FGMRES;

   // Preconditioner Options
   int prec_type = ifs_solver::PreconditionerType::SIMPLE; // enum PreconditionerType {SIMPLE = 0, SCR = 1, aSCR = 2, ILU = 3, Jacobi = 4};
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

   // track a node in the solid mesh
   bool track_node = false;
   FEParamVec3 m_track_node_pos;
   FILE *f = nullptr;
   Vector coord, disp;

   // output avg velocity of immersed solid
   bool output_avg_velocity = false;
   FILE *avg_velo_file = nullptr;
   Vector avg_velo;
   double volume = 0.0;

   // output additional immersed solid fields
   bool vonMisesStress = false;
public:
   // Custom stream buffer for writing to the log file
   CustomStreamBuffer customBuffer;

public:
   #ifdef MFEM_USE_MPI
   MFEMParFEMeshBuilder* par_femesh_builder = nullptr;
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
      // pre-partitioning of the mesh
      ADD_PARAMETER(m.pre_partition, "pre_partition_mesh");
      ADD_PARAMETER(m.partition_opt, "pre_partition_opt")->setEnums("PartitionMesh\0ReadParallelMesh\0");
      ADD_PARAMETER(m.np, FE_RANGE_GREATER_OR_EQUAL(1), "pre_partition_num_ranks");
      END_PARAM_GROUP();

   // Immersed Solid Options
   BEGIN_PARAM_GROUP("Immersed Solid Options");
      ADD_PARAMETER(m.immersed_solid, "Immersed_solid");
      ADD_PARAMETER(m.solidsolvertype, "Solid_solver_type")->setEnums("MFEM\0FEBio\0");
      ADD_PARAMETER(m.subtract_viscous_stress, "Subtract_viscous_stress");
      ADD_PARAMETER(m.three_field, "Three_field_formulation");
      ADD_PARAMETER(m.stabs_scale, FE_RANGE_GREATER_OR_EQUAL(1.0), "VMS_stabilization_scaling");
      ADD_PARAMETER(m.track_node, "Track_node");
      ADD_PARAMETER(m.m_track_node_pos, "Track_node_position");
      ADD_PARAMETER(m.output_avg_velocity, "Output_avg_velocity");
      ADD_PARAMETER(m.vonMisesStress, "Save_Von_Mises_Stress");
   END_PARAM_GROUP();

   // Time Stepping Options
   BEGIN_PARAM_GROUP("Time Stepping");
      ADD_PARAMETER(m.ts_method, "time_stepping_method")->setEnums("GenAlpha\0BackwardEuler\0");
      ADD_PARAMETER(m.rhoinf, FE_RANGE_CLOSED(0.0,1.0), "rho_inf");
      ADD_PARAMETER(m.step_save, FE_RANGE_GREATER_OR_EQUAL(1), "step_save");
      ADD_PARAMETER(m.steady_state_save, "steady_state_save");
      ADD_PARAMETER(m.steady_state_step_save, FE_RANGE_GREATER_OR_EQUAL(1), "steady_state_step_save");
      ADD_PARAMETER(m.steady_state_time, FE_RANGE_GREATER_OR_EQUAL(0.0), "steady_state_time");
      ADD_PARAMETER(m.write_restarts, "write_restart_files");
      ADD_PARAMETER(m.restart_step_save, FE_RANGE_GREATER_OR_EQUAL(1), "restart_step_save");
      ADD_PARAMETER(m.restarting, "restart_simulation");
      ADD_PARAMETER(m.restart_step, FE_RANGE_GREATER_OR_EQUAL(1), "restart_step");
      ADD_PARAMETER(m.write_final_state_as_IC, "write_final_state_as_IC");
      ADD_PARAMETER(m.final_state_folder, "final_state_folder");
      ADD_PARAMETER(m.start_from_specific_time, "start_from_specific_time");
      ADD_PARAMETER(m.start_sim_time, FE_RANGE_GREATER_OR_EQUAL(0.0), "start_sim_time");
      ADD_PARAMETER(m.start_sim_step, FE_RANGE_GREATER_OR_EQUAL(0.0), "start_sim_step");
      ADD_PARAMETER(m.wss, "Save_Wall_Shear_Stress");
      ADD_PARAMETER(m.vorticity, "Save_Vorticity");
      ADD_PARAMETER(m.q_criterion, "Save_Q_Criterion");
      ADD_PARAMETER(m.output_surface_pressures, "Output_surface_pressures");
      ADD_PARAMETER(m.output_surface_flowrates, "Output_surface_flow_rates");
   END_PARAM_GROUP();

   // // Nonlinear Solver Options
   BEGIN_PARAM_GROUP("Nonlinear Solver");
      ADD_PARAMETER(m.ns_no_error_nonconvergence, "no_error_nonconvergence");
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
      ADD_PARAMETER(m.prec_type, "prec_type")->setEnums("SIMPLE\0SCR\0aSCR\0ILU\0JACOBI\0");
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
   delete m.par_femesh_builder;
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

IncompressibleFlowSolver* MFEMSolver::GetFlowSolver()
{
   return m.flowsolver;
}

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

void MFEMSolver::AddWindkesselBC(FEParamDouble R, FEParamDouble Rd, FEParamDouble C, FEParamDouble p0, FEParamDouble pd, mfem::Array<int>& attr, bool compute_initial_pressure)
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
   double Q0 = m.flowsolver->ComputeMassFlowThruSurf(*m.flowsolver->GetCurrentVelocity(), attr, false, m.fluid_density, false);

   // compute initial pressure
   if (compute_initial_pressure)
   {
      double area;
      p0_val = m.flowsolver->ComputeAvgPressureOnSurf(*m.flowsolver->GetCurrentPressure(), attr, area);
   }   

   m.flowsolver->AddWindkesselBC(attr,R_val,Rd_val,C_val,p0_val,pd_val,Q0);
}

void MFEMSolver::AddMappedInflowBC_FromFile(std::string flowfile, bool periodic, int profile, mfem::Array<int>& attr)
{
   auto fluid_serial_mesh = m.immersed_solid ? static_cast<Mesh*>(m.fluid_mesh) : static_cast<Mesh*>(m.mesh);
   // TODO: for now just interp type 1 = SPLINE
   m.flowsolver->AddMappedInflowBC_FromFile(profile,fluid_serial_mesh,flowfile,false,attr,1,periodic);
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

   // Initialize the MFEM solver
   if (!CreateMesh()) return false;

   if (!SetupFluid()) return false;

   if(m.immersed_solid) if (!SetupSolid()) return false;

   SetupSolvers();

   PostProcessOutput();

   RestartFiles();

   #ifdef MFEM_USE_MPI
   if (!m.pre_partition)
   {
      ClearSerialMeshData();
   }
   #endif

   // solid post processing depends on parallel solid mesh
   if (m.immersed_solid) FinalizeSolid();

	return true;
}

bool MFEMSolver::CreateMesh()
{
   FEModel* fem = GetFEModel();

   // create the mesh
	FEMesh& femesh = fem->GetMesh();
	m.mesh = CreateMFEMMesh3D(femesh);
	if (m.mesh == nullptr) return false;
   
   #ifdef MFEM_USE_MPI
   if(m.pre_partition)
   {
      mfem::StopWatch timer_partitioning;
      timer_partitioning.Start();
      if (m.partition_opt == Imp::ParallelMesh::PartitionMesh)
      {
         // partition the mesh then output each piece
         if (m.immersed_solid)
         {
            // extract fluid and solid meshes
            if (!ExtractFluidAndSolidMeshes())
            {
               mfem::out << "Error extracting fluid and solid meshes" << std::endl;
               return false;
            }
         }

         mfem::out << "-----------------------------------------------------\n";
         mfem::out << "Writing mfem mesh partition data for " << std::to_string(m.np) << " ranks ..." << std::endl;
         mfem::out << "-----------------------------------------------------\n";

         // Partition the serial mesh and save the parallel mesh on np ranks
         // Create directory structure if it doesn't exist
         MFEM_ASSERT(mfem::Mpi::WorldSize() == 1, "Pre-partitioning should be done on a single rank (mpirun -np 1 ...)");
         #ifdef MFEM_USE_MPI
         int myid = mfem::Mpi::WorldRank();
         #else
         int myid = 0;
         #endif
         std::string part_dir = "ParallelMesh_np" + std::to_string(m.np) + "/";
         if (myid == 0)
         {
            try {
               std::filesystem::create_directories(part_dir);
            }
            catch (const std::filesystem::filesystem_error& ex) {
               MFEM_ABORT("Error creating directory for parallel mesh: " << ex.what());
            }
         
            std::string fluid_mname = part_dir + "fluid";
            Mesh* fluid_mesh_to_save = m.immersed_solid ? dynamic_cast<Mesh*>(m.fluid_mesh) : dynamic_cast<Mesh*>(m.mesh);
            SavePartitionedMesh(fluid_mname, *fluid_mesh_to_save, m.np);

            if (m.immersed_solid)
            {
               std::string solid_mname = part_dir + "solid";
               SavePartitionedMesh(solid_mname, *m.solid_mesh, m.np);
            }
         }        

         // This ensures that each BC has a MFEM attribute to write out
         PartitionMesh();
         if (SetupFluid() == false) return false;
         if (m.immersed_solid && m.solidsolvertype == SolidSolverType::FEBio)
         {
            if (SetupSolid() == false) return false;            
         }

         // write main and local .feb files
         std::string feb_filename = dynamic_cast<FEBioModel*>(fem)->GetFileTitle();
         feb_filename = std::filesystem::path(feb_filename).stem().string(); // remove .feb
         
         par_export.SetFEModel(fem);
         par_export.SetFESolver(this);
         bool immersed_febio_solid_and_solver = false;
         if (m.immersed_solid && m.solidsolvertype == SolidSolverType::FEBio)
         {
            immersed_febio_solid_and_solver = true;
            par_export.SetDomainAttributeMap(m.dom_attr_map);
         }
         
         mfem::out << "-----------------------------------------------------\n";
         mfem::out << "Writing parallel .feb files";
         par_export.WriteFiles(feb_filename, part_dir, m.np, immersed_febio_solid_and_solver);
         
         timer_partitioning.Stop();
         mfem::out << " ...SUCCESS!" << std::endl;
         mfem::out << "\nTime to partition: " << timer_partitioning.RealTime() << "s\n";
         mfem::out << "-----------------------------------------------------\n";
         std::exit(EXIT_SUCCESS);
      }
      else if (m.partition_opt == Imp::ParallelMesh::ReadParallelMesh)
      {
         mfem::out << "Reading the parallel mesh..." << std::endl;
         mfem::out << "-------------------------------\n";

         // read a parallel mesh from output above
         #ifdef MFEM_USE_MPI
         int myid = mfem::Mpi::WorldRank();
         #else
         int myid = 0;
         #endif
         std::string part_dir = "ParallelMesh_np" + std::to_string(mfem::Mpi::WorldSize()) + "/";

         std::string fluid_mname = part_dir + "fluid";
         m.fluid_pmesh = LoadParMesh(fluid_mname);

         // ensure attributes
         m.fluid_pmesh->SetAttributes();

         if (m.immersed_solid)
         {
            std::string solid_mname = part_dir + "solid";
            m.solid_pmesh = LoadParMesh(solid_mname);

            // ensure attributes
            m.solid_pmesh->SetAttributes();
         }

         // Read local feb file
         std::string feb_filename = dynamic_cast<FEBioModel*>(fem)->GetFileTitle();
         feb_filename = std::filesystem::path(feb_filename).stem().string(); // remove .feb
         std::string suffix = "_Main";
         if (feb_filename.size() >= suffix.size() &&
            feb_filename.compare(feb_filename.size() - suffix.size(), suffix.size(), suffix) == 0) {
            feb_filename.erase(feb_filename.size() - suffix.size());
         } // remove _Main
         std::ostringstream oss;
         oss << std::setw(10) << std::setfill('0') << myid;
         std::string local_feb_file = part_dir + feb_filename + "_Local_" + oss.str() + ".feb";

         if (m.immersed_solid) par_import.SetSolidMesh(m.solid_pmesh);
         if (par_import.Load(*fem, local_feb_file.c_str()) == false)
         {
            // Print the detailed import error captured by FEFileImport::errf
            char err[1024];
            par_import.GetErrorMessage(err);
            mfem::out << "MFEMParImport error (rank " << myid << "): "
                      << (err[0] ? err : "unknown error")
                      << "\nWhile loading: " << local_feb_file << std::endl;
            MFEM_ABORT("Error loading local .feb file");
         }

         // reinit model
         fem->Init();
      }
      else
      {
         MFEM_ABORT("Unknown option for handling parallel meshes");
      }
   }
   else
   {
      if (m.immersed_solid)
      {
         // extract fluid and solid meshes
         if (!ExtractFluidAndSolidMeshes())
         {
            mfem::out << "Error extracting fluid and solid meshes" << std::endl;
            return false;
         }
      }

      // Partition the meshes
      PartitionMesh();
   }   
   #endif

   // print mesh details
   if (m.immersed_solid)
   {
      feLog("-------------------------------\nFluid Mesh\n");
      #ifdef MFEM_USE_MPI
      m.fluid_pmesh->PrintInfo();
      #else
      m.fluid_mesh->PrintInfo();
      #endif
      feLog("-------------------------------\nSolid Mesh\n");
      #ifdef MFEM_USE_MPI
      m.solid_pmesh->PrintInfo();
      #else
      m.solid_mesh->PrintInfo();
      #endif
      feLog("-------------------------------\n");
   }
   else
   {
      feLog("-------------------------------\nFluid Mesh\n");
      #ifdef MFEM_USE_MPI
      m.fluid_pmesh->PrintInfo();
      #else
      m.mesh->PrintInfo();
      #endif
      feLog("-------------------------------\n");
   }

   return true;
}

bool MFEMSolver::ExtractFluidAndSolidMeshes()
{
   FEModel &fem = *GetFEModel();
   FEMesh &femesh = fem.GetMesh();

   // create fluid and solid submeshes from mfem version of FEBio mesh
   bool mfem_solid_material = false;
   bool febio_solid_material = false;
   Array<int> fluid_domain_attr(m.mesh->attributes.Max());
   fluid_domain_attr = 0;
   Array<int> solid_domain_attr(m.mesh->attributes.Max());
   solid_domain_attr = 0;
   for (int i = 0; i < femesh.Domains(); ++i)
   {
      FEDomain& dom = femesh.Domain(i);
      int dom_attr = m.mesh->GetDomainAttribute(&dom);
      if (dom.IsActive() && dynamic_cast<MFEMSolidMaterial*>(dom.GetMaterial()))
      {
         solid_domain_attr[i] = dom_attr;
         mfem_solid_material = true;
      }
      if (dom.IsActive() && dynamic_cast<FEElasticMaterial*>(dom.GetMaterial()))
      {
         solid_domain_attr[i] = dom_attr;
         febio_solid_material = true;

         FESolidDomain* solid_dom = dynamic_cast<FESolidDomain*>(&dom);
         MFEM_ASSERT(solid_dom != nullptr, "Expected solid domain");
         m.dom_attr_map[solid_dom] = dom_attr;
      }
      if(dom.IsActive() && dynamic_cast<MFEMNewtonianFluid*>(dom.GetMaterial()))
      {
         fluid_domain_attr[i] = dom_attr;            
      }
   }

   // create each mesh as sum of all domains (e.g. mfem solid mesh = sum of all febio solid domains), this will preserve domain numbering.
   // This assumes each solid domain is deformable.

   // fluid 
   m.fluid_mesh = new SubMesh(SubMesh::CreateFromDomain(*m.mesh,fluid_domain_attr));
   if (m.fluid_mesh == nullptr) return false;

   // solid
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

   return true;
}

#ifdef MFEM_USE_MPI
void MFEMSolver::SavePartitionedMesh(const std::string &prefix, Mesh &mesh, int np, int *partitioning)
{
   MeshPartitioner partitioner(mesh, np, partitioning);
   MeshPart mesh_part;
   for (int i = 0; i < np; i++)
   {
      partitioner.ExtractPart(i, mesh_part);

      ofstream f(MakeParFilename(prefix + "_", i, ".mesh", 10));
      f.precision(16);
      mesh_part.Print(f);
   }
}

ParMesh* MFEMSolver::LoadParMesh(const std::string &prefix)
{
   std::ifstream f(MakeParFilename(prefix + "_", Mpi::WorldRank(), ".mesh", 10));
   return new ParMesh(MPI_COMM_WORLD, f, false);
}

void MFEMSolver::PartitionMesh()
{
   auto fluid_serial_mesh = m.immersed_solid ? static_cast<Mesh*>(m.fluid_mesh) : static_cast<Mesh*>(m.mesh);
   m.fluid_pmesh = new ParMesh(MPI_COMM_WORLD, *fluid_serial_mesh);
   if (m.immersed_solid) m.solid_pmesh = new ParMesh(MPI_COMM_WORLD, *m.solid_mesh);

   if (m.immersed_solid && m.solidsolvertype == SolidSolverType::FEBio)
   {
      // create MFEMParFEMeshBuilder
      FEModel& fem = *GetFEModel();
      m.par_femesh_builder = new MFEMParFEMeshBuilder(&fem, m.mesh);
      m.par_femesh_builder->SetMFEMSubMesh(m.solid_mesh);
   }
}

void MFEMSolver::ClearSerialMeshData()
{
   FEModel* fem = GetFEModel();
	FEMesh& femesh = fem->GetMesh();

   // build local FEBio solid Mesh from mfem solid ParMesh
   femesh.Clear(); // clear FEBio serial mesh

   if (m.immersed_solid && m.solidsolvertype == SolidSolverType::FEBio)
   {
      ParFiniteElementSpace *dfes = m.flowsolver->GetImmersedSolidDisplacement()->ParFESpace();
      m.par_femesh_builder->BuildLocalFEBioMesh(dfes);

      // clear FEBio serial surface mesh data
      m.par_femesh_builder->ClearSerialSurfaces();
   }

   // delete the serial MFEM meshes
   delete m.mesh;
   delete m.fluid_mesh;
   delete m.solid_mesh;
}
#endif

bool MFEMSolver::SetupFluid()
{
   FEModel* fem = GetFEModel();
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
      m.flowsolver = new IncompressibleFlowSolver(m.fluid_pmesh, (int) m.orderv, (int) m.orderp, m.kinvis, m.fluid_density, m.orderingtype, m.immersed_solid);
      #else
      m.flowsolver = new IncompressibleFlowSolver(m.immersed_solid ? static_cast<Mesh*>(m.fluid_mesh) : static_cast<Mesh*>(m.mesh), (int) m.orderv, (int) m.orderp, m.kinvis, m.fluid_density, m.orderingtype, m.immersed_solid);
      #endif
   }
   
   // Time stepping method
   m.flowsolver->SetTimeSteppingMethod(m.ts_method);
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

   // flag for no error on non-convergence of newton solver
   // useful for immersed FSI with thin structures and vms stabilization 
   // scaling
   if (m.ns_no_error_nonconvergence) m.flowsolver->SetNewtonSolveNoErrorNonConvergence();

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
       case ifs_solver::StabilizationType::NONE:
            break;
		 case ifs_solver::StabilizationType::SUPG:
            if (m.orderv != m.orderp - 1) 
            {
                  feLogError("using SUPG only requires mixed order elements");
                  return false;
            }
            m.flowsolver->UseSUPG();
            break;

		 case ifs_solver::StabilizationType::SUPG_PSPG:
            if (m.orderv != m.orderp) 
            {
                  feLogError("SUPG/PSPG requires equal order elements");
                  return false;
            }
            m.flowsolver->UseSUPG();
            m.flowsolver->UsePSPG();
            break;

		 case ifs_solver::StabilizationType::VMS:
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
		m.flowsolver->SetViscousTermType(ifs_solver::ViscousTermType::LaplaceForm);
	}

   // Here we assume the immersed solid has no initial or boundary conditions set
   // Set the initial conditions.
	for (int i = 0; i < fem->InitialConditions(); ++i)
	{
		MFEMInitialCondition* ic = dynamic_cast<MFEMInitialCondition*>(fem->InitialCondition(i));
      if (dynamic_cast<MFEMInitialConditionFromFolder*>(ic)) continue; // skip ICs from folder, applied later on
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
         std::string surf_name;
         if (surf && surf->Elements() > 0)
         {
            FEFacetSet* set = surf->GetFacetSet();
            id = m.mesh->GetSurfaceAttribute(set); assert(id > 0);
            bc->SetBoundaryID(id);
            surf_name = surf->GetName();
            bc->SetSurfaceName(surf_name);
         }
         else
         {
            id = bc->GetBoundaryID(); assert(id > 0);
            surf_name = bc->GetSurfaceName();
         }		
         m.output_surfaces.push_back({surf_name, id});         
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
		if (sl && sl->IsActive() && !dynamic_cast<FEImmersedSurfaceLoad*>(sl))
		{
         FESurface* surf = &sl->GetSurface();
         int id = 0;
         std::string surf_name;
         if (surf && surf->Elements() > 0)
         {
            FEFacetSet* set = surf->GetFacetSet();
            id = m.mesh->GetSurfaceAttribute(set); assert(id > 0);
            sl->SetBoundaryID(id);
            surf_name = surf->GetName();
            sl->SetSurfaceName(surf_name);
         }
         else
         {
            id = sl->GetBoundaryID(); assert(id > 0);
            surf_name = sl->GetSurfaceName();
         }
         m.output_surfaces.push_back({surf_name, id}); 
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
         else if (bl->GetDomainAttr() != -1)
         {
            bl->Apply(this, bl->GetDomainAttr() - 1);
         }
         else
         {
            std::vector<int> dom_list;
            for (int j=0; j<bl->Domains(); ++j)
            {
               FEDomain* dom = bl->Domain(j);
               int id = m.mesh->GetDomainAttribute(dom); assert(id > 0);
               bl->Apply(this, id - 1);
               dom_list.push_back(id);
            }
            if (dom_list.size() == 1) bl->SetDomainAttr(dom_list[0]);
            else bl->SetDomainIDList(dom_list);
         }
      }
	}

   // Zero mean pressure
   if (m.zero_mean_pressure) m.flowsolver->AddZeroMeanPressureCondition();

   // Fluid grid functions
	m.u_gf = m.flowsolver->GetCurrentVelocity();
	m.p_gf = m.flowsolver->GetCurrentPressure();

   return true;
}

bool MFEMSolver::SetupSolid()
{
   FEModel *fem = GetFEModel();
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
      m.a_gf = m.flowsolver->GetImmersedSolidAcceleration();
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

      // add any immersed solid BCs
      for (int i = 0; i < fem->ModelLoads(); ++i)
      {
         // surface loads
         FEImmersedSurfaceLoad* im_sl = dynamic_cast<FEImmersedSurfaceLoad*>(fem->ModelLoad(i));
         if (im_sl && im_sl->IsActive())
         {
            FESurface* surf = &im_sl->GetSurface();
            int id = 0;
            if (surf)
            {
               FEFacetSet* set = surf->GetFacetSet();
               id = m.mesh->GetSurfaceAttribute(set); assert(id > 0);
               im_sl->SetBoundaryID(id); // store for later use
            }
            else
            {
               id = im_sl->GetBoundaryID(); assert(id > 0);
            }

            im_sl->Apply(this, id - 1);
            m.immersed_febio_solid->immersed_surface_loads = true;
         }
      }

      // immersed solid grid functions
      m.d_gf = m.flowsolver->GetImmersedSolidDisplacement();
      m.vib_gf = m.flowsolver->GetImmersedSolidVelocity();
      m.a_gf = m.flowsolver->GetImmersedSolidAcceleration();
   }
   else
   {
      feLogError("Invalid immersed solid solver type");
      return false;
   }      

   // Create transfer maps and global functions 
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

   return true;
}

void MFEMSolver::SetupSolvers()
{
   FEModel *fem = GetFEModel();
   double dt = fem->GetTime().timeIncrement;
	m.flowsolver->Setup(dt);
}

void MFEMSolver::PostProcessOutput()
{
   // setup output in paraview format
   m.fname = "results_mfem_solver";
   int curr_step;
   double curr_time;
   if (m.restarting)
   {
      m.fname += "_RESTART_STEP_" + std::to_string(m.restart_step);
      curr_step = m.restart_step;
      curr_time = m.restart_t;
   }
   else if (m.start_from_specific_time)
   {
      m.fname += "_START_STEP_" + std::to_string(m.start_sim_step);
      curr_step = m.start_sim_step;
      curr_time = m.start_sim_time;
   }
   else
   {
      curr_step = 0;
      curr_time = 0.0;
   }
   #ifdef MFEM_USE_MPI
   m.fname += "_np" + std::to_string(mfem::Mpi::WorldSize());
   #endif

   std::string fluid_filename = m.fname + "/fluid";
   #ifdef MFEM_USE_MPI
   #ifdef MFEM_USE_HDF5
   m.pvdc = new ParaViewHDFDataCollection(fluid_filename, m.fluid_pmesh);
   #else
   m.pvdc = new ParaViewDataCollection(fluid_filename, m.fluid_pmesh);
   #endif
   #else
   m.pvdc = new ParaViewDataCollection(fluid_filename, m.fluid_mesh);
   #endif   
	m.pvdc->SetDataFormat(m.output_format);
	m.pvdc->SetHighOrderOutput(false);
	m.pvdc->SetLevelsOfDetail(m.orderv);
   m.pvdc->SetCompression(true);
   m.pvdc->SetCompressionLevel(m.compression_level);
	m.pvdc->SetCycle(curr_step);
	m.pvdc->SetTime(curr_time);
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
   if (m.wss)
   {
      m.flowsolver->ComputeWSS(*m.u_gf,*m.p_gf);
      m.wss_gf = m.flowsolver->GetCurrentWSS();
      m.pvdc->RegisterField("WSS", m.wss_gf);
   }
   if (m.vorticity)
   {
      m.vorticity_gf.SetSpace(m.u_gf->FESpace());
      m.vorticity_gf = 0.0;
      m.flowsolver->ComputeCurl3D(*m.u_gf,m.vorticity_gf);
      m.pvdc->RegisterField("Vorticity", &m.vorticity_gf);
   }
   if (m.q_criterion)
   {
      m.q_criterion_gf.SetSpace(m.u_gf->FESpace());
      m.q_criterion_gf = 0.0;
      ComputeQCriterion(*m.u_gf,m.q_criterion_gf);
      m.pvdc->RegisterField("Q-Criterion", &m.q_criterion_gf);
   }
	m.pvdc->Save();

   if (m.immersed_solid)
   {
      std::string solid_filename = m.fname + "/solid";
      #ifdef MFEM_USE_MPI
      #ifdef MFEM_USE_HDF5
      m.sdc = new ParaViewHDFDataCollection(solid_filename, m.solid_pmesh);
      #else
      m.sdc = new ParaViewDataCollection(solid_filename, m.solid_pmesh);
      #endif
      #else
      m.sdc = new ParaViewDataCollection(solid_filename, m.solid_mesh);  
      #endif
      m.sdc->SetDataFormat(m.output_format);
      m.sdc->SetHighOrderOutput(false);
      m.sdc->SetLevelsOfDetail(m.orderd);
      m.sdc->SetCompression(true);
      m.sdc->SetCompressionLevel(m.compression_level);
      m.sdc->RegisterField("Displacement", m.d_gf);
      m.sdc->RegisterField("Velocity", m.vib_gf);
      m.sdc->RegisterField("Acceleration", m.a_gf);

      m.sdc->SetCycle(curr_step);
      m.sdc->SetTime(curr_time);
   }

   if (m.output_surface_pressures)
   {
      #ifdef MFEM_USE_MPI
      m.surf_press_ts = std::make_unique<SurfacePressureTimeSeries>(
         m.fluid_pmesh, m.p_gf, m.flowsolver);
      #else
      m.surf_press_ts = std::make_unique<SurfacePressureTimeSeries>(
         m.fluid_mesh, m.p_gf, m.flowsolver);
      #endif

      FEModel* fem = GetFEModel();
      FEMesh &mesh = fem->GetMesh();
      if (mesh.Surfaces() == 0)
      {
         for (int i=0; i<m.output_surfaces.size(); ++i) 
         {
            std::string sname = m.output_surfaces[i].name;
            int attr = m.output_surfaces[i].attribute;
            m.surf_press_ts->AddSurface(sname, attr);
         }
      }
      else
      {
         for (int i = 0; i < mesh.Surfaces(); ++i)
         {
            FESurface& surf = mesh.Surface(i);
            if (surf.IsActive())
            {
               int attr = m.mesh->GetSurfaceAttribute(surf.GetFacetSet());
               std::string sname = surf.GetName().empty()
                                 ? ("Surface_" + std::to_string(attr))
                                 : surf.GetName();
               m.surf_press_ts->AddSurface(sname, attr);
            }
         }
      }
      if (m.surf_press_ts->Initialize(m.fname, "surface_pressures.txt"))
      {
         double start_time = m.restarting ? m.restart_t :
                           (m.start_from_specific_time ? m.start_sim_time : 0.0);
         m.surf_press_ts->Append(start_time);
      }
   }
   if (m.output_surface_flowrates) 
   {
      #ifdef MFEM_USE_MPI
      m.surf_flow_ts = std::make_unique<SurfaceFlowRateTimeSeries>(
         m.fluid_pmesh, m.u_gf, m.flowsolver, m.fluid_density);
      #else
      m.surf_flow_ts = std::make_unique<SurfaceFlowRateTimeSeries>(
         m.fluid_mesh, m.u_gf, m.flowsolver, m.fluid_density);
      #endif

      FEModel* fem = GetFEModel();
      FEMesh &mesh = fem->GetMesh();
      if (mesh.Surfaces() == 0) 
      {
         for (int i=0; i<m.output_surfaces.size(); ++i) 
         {
            std::string sname = m.output_surfaces[i].name;
            int attr = m.output_surfaces[i].attribute;
            m.surf_flow_ts->AddSurface(sname, attr);
         }
      }
      else 
      {
         for (int i = 0; i < mesh.Surfaces(); ++i) 
         {
            FESurface& surf = mesh.Surface(i);
            if (surf.IsActive()) 
            {
               int attr = m.mesh->GetSurfaceAttribute(surf.GetFacetSet());
               std::string sname = surf.GetName().empty()
                                    ? ("Surface_" + std::to_string(attr))
                                    : surf.GetName();
               m.surf_flow_ts->AddSurface(sname, attr);
            }
         }
      }
      if (m.surf_flow_ts->Initialize(m.fname, "surface_flow_rates.txt")) 
      {
         double start_time = m.restarting ? m.restart_t :
                           (m.start_from_specific_time ? m.start_sim_time : 0.0);
         m.surf_flow_ts->Append(start_time);
      }
   }
}

void MFEMSolver::RestartFiles()
{
   // Setup restart files
   if (m.write_restarts)
   {
      m.flowsolver->writeRestartFiles(m.restart_folder, m.restart_step_save);
   }
   if (m.restarting)
   {
      m.flowsolver->Restarting(m.restart_folder, m.restart_step, m.restart_dt, m.restart_t);
      m.flowsolver->readRestarts(m.restart_step, m.restart_dt, m.restart_t);

      // save the restart state
      if (m.wss) m.flowsolver->ComputeWSS(*m.u_gf,*m.p_gf);
      if (m.vorticity) m.flowsolver->ComputeCurl3D(*m.u_gf, m.vorticity_gf);
      if (m.q_criterion) ComputeQCriterion(*m.u_gf, m.q_criterion_gf);

      m.pvdc->SetCycle(m.restart_step);
      m.pvdc->SetTime(m.restart_t);
      m.pvdc->Save();
      if (m.immersed_solid)
      {
         m.sdc->SetCycle(m.restart_step);
         m.sdc->SetTime(m.restart_t);
      }
   }
   if (m.start_from_specific_time)
   {
      MFEM_ASSERT(!m.restarting, "Cannot restart and start from specific time");
      FEModel *fem = GetFEModel();

      for (int i = 0; i < fem->InitialConditions(); ++i)
      {
         MFEMInitialCondition* ic = dynamic_cast<MFEMInitialCondition*>(fem->InitialCondition(i));
         if (dynamic_cast<MFEMInitialConditionFromFolder*>(ic) && ic->IsActive())
         {
            ic->Apply(this);
         }
      }

      // save the state
      if (m.wss) m.flowsolver->ComputeWSS(*m.u_gf,*m.p_gf);
      if (m.vorticity) m.flowsolver->ComputeCurl3D(*m.u_gf, m.vorticity_gf);
      if (m.q_criterion) ComputeQCriterion(*m.u_gf, m.q_criterion_gf);

      m.pvdc->SetCycle(m.start_sim_step);
      m.pvdc->SetTime(m.start_sim_time);
      m.pvdc->Save();
      if (m.immersed_solid)
      {
         m.sdc->SetCycle(m.start_sim_step);
         m.sdc->SetTime(m.start_sim_time);
      }
   }
}

void MFEMSolver::FinalizeSolid()
{
   // ensure new febio solid domain has proper data after building parallel mesh
   if (m.immersed_solid && m.solidsolvertype == SolidSolverType::FEBio) m.immersed_febio_solid->SetupFEDomains();

   // additional solid post processing output
   if (m.vonMisesStress)
   {
      auto immersed_solid_solver = m.flowsolver->GetImmersedSolid();
      m.vonMisesStress_gf = immersed_solid_solver->ComputeNodallyAveragedVonMisesStress(*m.d_gf);
      m.sdc->RegisterField("VonMisesStress", m.vonMisesStress_gf);
   }

   m.sdc->Save();
}

bool MFEMSolver::SolveStep()
{
	// Get the FEBio model and mesh objects
	FEModel& fem = *GetFEModel();
	FEMesh& femesh = fem.GetMesh();

	FEAnalysis* analysis = fem.GetCurrentStep();
   double t, dt;
   int timeStep;
   if (m.restarting)
   {
      dt = m.restart_dt;
      t = m.restart_t + dt; // t = t_n+1
      timeStep = m.restart_step + 1;

      // update FEBio time step data
      analysis->m_ntimesteps = m.restart_step; // time steps completed
      analysis->m_dt = dt;

      fem.GetTime().currentTime = t;
      fem.GetTime().timeIncrement = dt;
      fem.GetTime().timeStep = dt;

      m.restarting = false; // only restart once
   }
   else if (m.start_from_specific_time)
   {
      dt = fem.GetTime().timeIncrement;

      // check that start time and start step are compatible with dt
      double start_time = dt * (double) m.start_sim_step;
      MFEM_ASSERT(start_time == m.start_sim_time, "Error: start time and start step are not compatible with dt");

      t = m.start_sim_time + dt; // t = t_n+1
      timeStep = m.start_sim_step + 1; // current time step number

      // update FEBio time step data
      analysis->m_ntimesteps = m.start_sim_step; // time steps completed
      analysis->m_dt = dt;

      fem.GetTime().currentTime = t;
      fem.GetTime().timeIncrement = dt;
      fem.GetTime().timeStep = dt;

      m.start_from_specific_time = false; // only start from specific time once
   }
   else
   {
      t = fem.GetTime().currentTime; // t = t_n+1
      dt = fem.GetTime().timeIncrement;
      timeStep = analysis->m_ntimesteps + 1; // current time step number
   }

   // solve a time step
	m.flowsolver->Step(t, dt, timeStep);
   
	// save data
   int curr_step_save = m.step_save;
   if (m.steady_state_save && t >= m.steady_state_time)
   {
      curr_step_save = m.steady_state_step_save;
   }   
   if ((timeStep) % curr_step_save == 0)
   {
      m.flowsolver->sw_io.Start();
      if (m.wss) m.flowsolver->ComputeWSS(*m.u_gf,*m.p_gf);
      if (m.vorticity) m.flowsolver->ComputeCurl3D(*m.u_gf, m.vorticity_gf);
      if (m.q_criterion) ComputeQCriterion(*m.u_gf, m.q_criterion_gf);

      m.pvdc->SetCycle(timeStep);
      m.pvdc->SetTime(t);
      m.pvdc->Save();
      if (m.immersed_solid)
      {
         if (m.vonMisesStress)
         {
            auto immersed_solid_solver = m.flowsolver->GetImmersedSolid();
            m.vonMisesStress_gf = immersed_solid_solver->ComputeNodallyAveragedVonMisesStress(*m.d_gf);
         }

         m.sdc->SetCycle(timeStep);
         m.sdc->SetTime(t);
         m.sdc->Save();
      }
      m.flowsolver->sw_io.Stop();
   }	

   if (m.output_surface_pressures && m.surf_press_ts)
   {
      m.flowsolver->sw_io.Start();
      m.surf_press_ts->Append(t);
      m.flowsolver->sw_io.Stop();
   }
   if (m.output_surface_flowrates && m.surf_flow_ts) 
   {
      m.flowsolver->sw_io.Start();
      m.surf_flow_ts->Append(t);
      m.flowsolver->sw_io.Stop();
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

   // store data in txt files
   if (m.immersed_solid && m.track_node)
   {
      m.flowsolver->sw_io.Start();
      auto immersed_solidsolver = m.flowsolver->GetImmersedSolid();
      immersed_solidsolver->TrackNodalDisplacement(m.coord,m.disp);

      #ifdef MFEM_USE_MPI
      if(Mpi::Root()) {
      #endif
         fprintf(m.f, "%20.16e     %20.16e      %20.16e      %20.16e     %20.16e\n",t,m.disp(0),m.disp(1),m.disp(2),m.disp.Norml2());
         fflush(m.f);
         fflush(stdout);
      #ifdef MFEM_USE_MPI
      }
      #endif
      m.flowsolver->sw_io.Stop();
   }
   if (m.immersed_solid && m.output_avg_velocity)
   {
      m.flowsolver->sw_io.Start();
      ComputeAverageVelocity(m.form_type,*m.d_gf,*m.vib_gf,m.avg_velo,m.volume);

      #ifdef MFEM_USE_MPI
      if(Mpi::Root()) {
      #endif
         fprintf(m.avg_velo_file, "%20.16e     %20.16e     %20.16e     %20.16e     %20.16e\n",t,m.avg_velo(0),m.avg_velo(1),m.avg_velo(2),m.volume);  
         fflush(m.avg_velo_file);
         fflush(stdout);
      #ifdef MFEM_USE_MPI
      }
      #endif
      m.flowsolver->sw_io.Stop();
   }

	// print out mfem timing data
	if (timeStep == analysis->m_ntime)
	{
      // timing info
		m.flowsolver->PrintTimingData();

      // write final state as restart files for future use
      if (m.write_final_state_as_IC)
      {
         m.flowsolver->writeRestartFiles(m.final_state_folder,1); // second arg not used here
         m.flowsolver->writeRestarts(0,0,0);
      }

      // // print out relative L2 error at final time step
      // double err_u = 0.0;
      // double err_p = 0.0;

      // // exact solutions
      // auto vel_exact = [](const Vector &x, double t, Vector &u) -> void {
      //    u.SetSize(x.Size());
      //    u = 0.0;
      // };

      // auto grad_vel_exact = [](const Vector &x, double t, DenseMatrix &gradu) -> void {
      //    gradu.SetSize(x.Size());
      //    gradu = 0.0;
      // };

      // auto pres_exact = [](const Vector &x, double t) -> double {
      //    double p_return = 0.0;

      //    // parameters
      //    double R = 0.25;
      //    double w = 0.06250;
      //    double mu = 1.0;
      //    double len = 1.0;

      //    // compute r
      //    Vector x_shifted(x.Size());
      //    x_shifted = x;
      //    x_shifted(0) -= 0.5;
      //    x_shifted(1) -= 0.5;

      //    // compute magnitude of vector
      //    double r = x_shifted.Norml2();

      //    // compute pres
      //    if (r >= R + w)
      //    {
      //       p_return = - M_PI * mu * ( (R + w) * (R + w) - R * R)/(2.0 * len*len);
      //    }
      //    else if (r > R && r < (R + w))
      //    {
      //       p_return = mu * log((R + w)/r) - M_PI * mu * ( (R + w) * (R + w) - R * R)/(2.0 * len*len);
      //    }
      //    else if (r <= R)
      //    {
      //       p_return = mu * log(1.0 + w/R) - M_PI * mu * ( (R + w) * (R + w) - R * R)/(2.0 * len*len);
      //    }
      //    else
      //    {
      //       MFEM_ABORT("Error: r is out of range");
      //    }

      //    return p_return;
      // };

      // VectorFunctionCoefficient vel_excoeff(m.flowsolver->GetMesh()->Dimension(), vel_exact);
      // FunctionCoefficient p_excoeff(pres_exact);
      // vel_excoeff.SetTime(t);
      // p_excoeff.SetTime(t);
      
      // // compute relative L2 error on finer mesh
      // int order_err = m.orderv + 3;
      // const IntegrationRule *irs[Geometry::NumGeom];
      // for (int i=0; i < Geometry::NumGeom; ++i)
      // {
      //    irs[i] = &(IntRules.Get(i, order_err));
      // }
      // err_u  = m.u_gf->ComputeL2Error(vel_excoeff, irs);
      // err_p  = m.p_gf->ComputeL2Error(p_excoeff, irs);

      // // compute velocity h1 norm
      // MatrixFunctionCoefficient gradu_excoeff(m.flowsolver->GetMesh()->Dimension(),grad_vel_exact);

      // // set time on coefficient
      // gradu_excoeff.SetTime(t);
      // double h1_err_u = ComputeVectorH1Error(m.u_gf,&vel_excoeff,&gradu_excoeff,irs);

      // #ifdef MFEM_USE_MPI
      // if (Mpi::Root())
      // {
      //    mfem::out << setprecision(16);
      //    // absolute h1 error
      //    // h1 norm for velocity
      //    mfem::out << "|| u_h - u_ex ||_h1 = " << h1_err_u << "\n";

      //    // absolute L2 error
      //    mfem::out << "|| u_h - u_ex || = " << err_u << "\n";
      //    mfem::out << "|| p_h - p_ex || = " << err_p << "\n";  
      // }
      // #else
      // mfem::out << setprecision(16);
      // // absolute h1 error
      // // h1 norm for velocity
      // mfem::out << "|| u_h - u_ex ||_h1 = " << h1_err_u << "\n";

      // // absolute L2 error
      // mfem::out << "|| u_h - u_ex || = " << err_u << "\n";
      // mfem::out << "|| p_h - p_ex || = " << err_p << "\n";  
      // #endif
	}

	return true;
}

double ComputeVectorH1Error(mfem::GridFunction *u_gf, mfem::VectorCoefficient *ex_u, mfem::MatrixCoefficient *ex_gradu, const mfem::IntegrationRule *irs[])
{
   double error_temp = 0.0;
   const mfem::FiniteElement *fe;
   mfem::ElementTransformation *Tr;
   mfem::DenseMatrix grad;
   int intorder;
   int dim = u_gf->FESpace()->GetMesh()->SpaceDimension();
   mfem::DenseMatrix mat(dim);
   Vector uvec(dim), vec(dim);

   for (int i = 0; i < u_gf->FESpace()->GetNE(); i++)
   {
      fe = u_gf->FESpace()->GetFE(i);
      Tr = u_gf->FESpace()->GetElementTransformation(i);
      intorder = 2*fe->GetOrder() + 3; // <--------
      const mfem::IntegrationRule *ir;
      if (irs)
      {
         ir = irs[fe->GetGeomType()];
      }
      else
      {
         ir = &(IntRules.Get(fe->GetGeomType(), intorder));
      }
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(j);
         Tr->SetIntPoint(&ip);
         u_gf->GetVectorGradient(*Tr,grad);
         ex_gradu->Eval(mat,*Tr,ip);
         mat.Add(-1.0,grad);

         u_gf->GetVectorValue(*Tr,ip,vec);
         ex_u->Eval(uvec,*Tr,ip);
         vec -= uvec;

         error_temp += ip.weight * Tr->Weight() * ( mat.FNorm2() + vec*vec);
      }
   }

   #ifdef MFEM_USE_MPI
   // parallelization
   double glob_error = 0.0;
   double local_error = abs(error_temp);
   MPI_Allreduce(&local_error, &glob_error, 1, MPI_DOUBLE, MPI_SUM, dynamic_cast<mfem::ParGridFunction*>(u_gf)->ParFESpace()->GetComm());
   #else
   double glob_error = abs(error_temp);
   #endif
   
   return sqrt(glob_error);
}

// Computes Q = 0.5*(tr(\nabla u)^2 - tr(\nabla u \cdot \nabla u))
void ComputeQCriterion(GridFunction &u_gf, GridFunction &q_gf)
{
   FiniteElementSpace *v_fes = u_gf.FESpace();
   FiniteElementSpace *fes = q_gf.FESpace();

   // AccumulateAndCountZones
   Array<int> zones_per_vdof;
   zones_per_vdof.SetSize(fes->GetVSize());
   zones_per_vdof = 0;

   q_gf = 0.0;

   // Local interpolation
   int elndofs;
   Array<int> v_dofs, dofs;
   Vector vals;
   Vector loc_data;
   int vdim = v_fes->GetVDim();
   DenseMatrix grad_hat;
   DenseMatrix dshape;
   DenseMatrix grad;

   for (int e = 0; e < fes->GetNE(); ++e)
   {
      fes->GetElementVDofs(e, dofs);
      v_fes->GetElementVDofs(e, v_dofs);
      u_gf.GetSubVector(v_dofs, loc_data);
      vals.SetSize(dofs.Size());
      ElementTransformation *tr = fes->GetElementTransformation(e);
      const FiniteElement *el = fes->GetFE(e);
      elndofs = el->GetDof();
      int dim = el->GetDim();
      dshape.SetSize(elndofs, dim);

      for (int dof = 0; dof < elndofs; ++dof)
      {
         // Project
         const IntegrationPoint &ip = el->GetNodes().IntPoint(dof);
         tr->SetIntPoint(&ip);

         // Eval
         // GetVectorGradientHat
         el->CalcDShape(tr->GetIntPoint(), dshape);
         grad_hat.SetSize(vdim, dim);
         DenseMatrix loc_data_mat(loc_data.GetData(), elndofs, vdim);
         MultAtB(loc_data_mat, dshape, grad_hat);

         const DenseMatrix &Jinv = tr->InverseJacobian();
         grad.SetSize(grad_hat.Height(), Jinv.Width());
         Mult(grad_hat, Jinv, grad);

         double q_val = 0.5 * (sq(grad(0, 0)) + sq(grad(1, 1)) + sq(grad(2, 2)))
                        + grad(0, 1) * grad(1, 0) + grad(0, 2) * grad(2, 0)
                        + grad(1, 2) * grad(2, 1);

         vals(dof) = q_val;
      }

      // Accumulate values in all dofs, count the zones.
      for (int j = 0; j < dofs.Size(); j++)
      {
         int ldof = dofs[j];
         q_gf(ldof) += vals[j];
         zones_per_vdof[ldof]++;
      }
   }

   #ifdef MFEM_USE_MPI
   // Communication

   // Count the zones globally.
   mfem::GroupCommunicator &gcomm = dynamic_cast<ParGridFunction*>(&q_gf)->ParFESpace()->GroupComm();
   gcomm.Reduce<int>(zones_per_vdof, mfem::GroupCommunicator::Sum);
   gcomm.Bcast(zones_per_vdof);

   // Accumulate for all vdofs.
   gcomm.Reduce<double>(q_gf.GetData(), mfem::GroupCommunicator::Sum);
   gcomm.Bcast<double>(q_gf.GetData());

   // Compute means
   for (int i = 0; i < q_gf.Size(); i++)
   {
      const int nz = zones_per_vdof[i];
      if (nz)
      {
         q_gf(i) /= nz;
      }
   }
   #endif
}

void ComputeAverageVelocity(int form_type, GridFunction &u_gf, GridFunction &v_gf, Vector &avg_velo, double &volume)
{
   const FiniteElement *fe;
   ElementTransformation *Tr;
   DenseMatrix gradu,F,I;
   int intorder;
   int dim = v_gf.FESpace()->GetMesh()->SpaceDimension();   
   gradu.SetSize(dim);
   gradu = 0.0;
   F.SetSize(dim);
   F = 0.0;
   I.Diag(1.0,dim);
   double J = 0.0;
   double vx = 0.0;
   double vy = 0.0;
   double vz = 0.0;
   double volume_integral = 0.0;
   Vector v(dim);

   for (int i = 0; i < v_gf.FESpace()->GetNE(); i++)
   {
      fe = v_gf.FESpace()->GetFE(i);
      Tr = v_gf.FESpace()->GetElementTransformation(i);
      intorder = 2*fe->GetOrder() + 3; // <--------
      const IntegrationRule *ir;
      ir = &(IntRules.Get(fe->GetGeomType(), intorder));
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir->IntPoint(j);
         Tr->SetIntPoint(&ip);

         if (form_type == ImmersedSolid::form_type::TotalLagrangian)
         {
            // Compute J
            u_gf.GetVectorGradient(*Tr,gradu);
            Add(1.0,gradu,1.0,I,F);
            J = F.Det();

            // Compute v
            v_gf.GetVectorValue(*Tr,ip,v);

            // v integral
            vx += ip.weight * Tr->Weight() * v(0) * J;
            vy += ip.weight * Tr->Weight() * v(1) * J;
            if (dim == 3) vz += ip.weight * Tr->Weight() * v(2) * J;

            // volume integral
            volume_integral += ip.weight * Tr->Weight() * J;
         }
         else
         {
            // Compute v
            v_gf.GetVectorValue(*Tr,ip,v);

            // v integral
            vx += ip.weight * Tr->Weight() * v(0);
            vy += ip.weight * Tr->Weight() * v(1);
            if (dim == 3) vz += ip.weight * Tr->Weight() * v(2);

            // volume integral
            volume_integral += ip.weight * Tr->Weight();
         }      
      }
   }

   #ifdef MFEM_USE_MPI
   // parallelization
   double glob_vx = 0.0;
   double glob_vy = 0.0;
   double glob_vz = 0.0;
   double glob_volume = 0.0;
   MPI_Allreduce(&vx, &glob_vx, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, dynamic_cast<ParGridFunction*>(&u_gf)->ParFESpace()->GetComm());
   MPI_Allreduce(&vy, &glob_vy, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, dynamic_cast<ParGridFunction*>(&u_gf)->ParFESpace()->GetComm());
   if (dim == 3) MPI_Allreduce(&vz, &glob_vz, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, dynamic_cast<ParGridFunction*>(&u_gf)->ParFESpace()->GetComm());
   MPI_Allreduce(&volume_integral, &glob_volume, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, dynamic_cast<ParGridFunction*>(&u_gf)->ParFESpace()->GetComm());
   #else
   double glob_vx = vx;
   double glob_vy = vy;
   double glob_vz = vz;
   double glob_volume = volume_integral;
   #endif

   avg_velo.SetSize(dim);
   avg_velo(0) = glob_vx / glob_volume;
   avg_velo(1) = glob_vy / glob_volume;
   if (dim == 3) avg_velo(2) = glob_vz / glob_volume;

   volume = glob_volume;
}