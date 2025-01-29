#include "MFEMSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#ifdef WIN32
#include "ins_solver.hpp"
#include "ins_solver_enum.hpp"
#include "immersed_solid.hpp" // file in ImmersedFSI-Serial/
#include "materials.hpp" // file in Elasticity/Materials/
#else
#include <Navier-Stokes-Serial/ins_solver.hpp>
#include <Navier-Stokes/ins_solver_enum.hpp>
#include <ImmersedFSI-Serial/immersed_solid.hpp>
#include <Elasticity/Materials/materials.hpp>
#include "mfem.hpp"
#endif
#include "MFEMMesh.h"
#include "MFEMInitialCondition.h"
#include "MFEMBoundaryCondition.h"
#include "MFEMLoad.h"
#include <FECore/log.h>
#include "MFEMMaterial.h"
#include "CustomStreamBuffer.h"
#include <typeinfo>

using namespace ins_solver_serial;
using namespace std;
using namespace std::placeholders;
using namespace mfem;
using namespace immersed_solid_serial;
using namespace mat;

class MFEMSolver::Imp
{
public:
	static Imp* This;
	Imp() { This = this; }

public:
	MFEMMesh* mesh = nullptr; // MFEM mesh
	INS_Solver* flowsolver = nullptr; // MFEM flow solver
	GridFunction* p_gf = nullptr; // velocity grid function
	GridFunction* u_gf = nullptr; // pressure grid function
	ParaViewDataCollection* pvdc = nullptr; // paraview data collection

   MFEMMesh* solid_mesh = nullptr; // immersed solid mesh
   GridFunction *d_gf = nullptr; // immersed solid displacement grid function
   GridFunction *vib_gf = nullptr; // immersed solid velocity grid function
   ParaViewDataCollection *sdc = nullptr; // immersed solid paraview data collection

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

   int nonlinear_form = ins_solver::NonlinearFormType::Convective;
   int viscous_term_type = ins_solver::ViscousTermType::StressForm;
   bool LaplaceForm = false;
   bool total_integ = true;
   bool zero_mean_pressure = false;

   // Time Stepping Options
   double rhoinf = 0.5;

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

   // Preconditioner Options
   int prec_type = ins_solver::PreconditionerType::SIMPLE;
   bool direct_solver_prec = false;
   double prec_Ainv_atol = 1.0e-6;
   double prec_Ainv_rtol = 1.0e-3;
   double prec_Sinv_atol = 1.0e-6;
   double prec_Sinv_rtol = 1.0e-3;

   // immersed solid
   int orderd = orderv;
   MaterialType *mat_type = nullptr;
   
   int form_type = ImmersedSolid::form_type::TotalLagrangian;
   // int form_type = ImmersedSolid::form_type::UpdatedLagrangian; // either work

   bool subtract_viscous_stress = true;
   bool immersed_solid = false;
   int vol_type = NearlyIncompressibleMaterial::volumetric_energy::ST91;
   // TO DO make vol_type a drop down option with nearly incompressible materials
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
      ADD_PARAMETER(m.subtract_viscous_stress, "Subtract_viscous_stress");
      // drop down menu to select Total vs Updated Lagrangian TODO
   END_PARAM_GROUP();

   // Time Stepping Options
   BEGIN_PARAM_GROUP("Time Stepping");
      ADD_PARAMETER(m.rhoinf, FE_RANGE_CLOSED(0.0,1.0), "rho_inf");
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
      ADD_PARAMETER(m.direct_solver, "direct_solver");
      ADD_PARAMETER(m.lin_atol, FE_RANGE_GREATER_OR_EQUAL(0.0), "lin_atol");
      ADD_PARAMETER(m.lin_rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "lin_rtol");
      // ADD_PARAMETER(m.lin_minIt, FE_RANGE_GREATER_OR_EQUAL(1), "lin_minIt");
      ADD_PARAMETER(m.lin_maxIt, FE_RANGE_GREATER_OR_EQUAL(1), "lin_maxIt");
   END_PARAM_GROUP();

   // Preconditioner Options
   BEGIN_PARAM_GROUP("Preconditioner");
      // set preconditioner drop down menu, for now SIMPLE
      ADD_PARAMETER(m.direct_solver_prec, "inner_solves_using_direct_solver");
      ADD_PARAMETER(m.prec_Ainv_atol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ainv_atol");
      ADD_PARAMETER(m.prec_Ainv_rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ainv_rtol");
      ADD_PARAMETER(m.prec_Sinv_atol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Sinv_atol");
      ADD_PARAMETER(m.prec_Sinv_rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Sinv_rtol");
   END_PARAM_GROUP();
END_FECORE_CLASS();

MFEMSolver::MFEMSolver(FEModel* fem) : FESolver(fem), m(*(new MFEMSolver::Imp))
{
   // default values created in Imp class
}

MFEMSolver::~MFEMSolver()
{
   // Restore the original buffer
   mfem::out.rdbuf(customBuffer.originalBuffer);

	delete m.mesh;
	delete m.flowsolver;
	delete m.pvdc;

   delete m.solid_mesh;
   delete m.mat_type;
   delete m.sdc;
}

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
   GridFunction *u_ic = m.flowsolver->GetCurrentVelocity();
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
	FEModel& fem = *GetFEModel();

   // Redirect mfem::out to the custom stream buffer that writes to the log file
   // TODO: Cleanup the output to the log file
   customBuffer.originalBuffer = mfem::out.rdbuf();
   customBuffer.fe_model = &fem;
   std::ostream customStream(&customBuffer);
   mfem::out.rdbuf(customStream.rdbuf());

   // create fluid mesh
	FEMesh& femesh = fem.GetMesh();
	m.mesh = CreateMFEMMesh3D(femesh);   
	if (m.mesh == nullptr) return false;
   if (m.immersed_solid)
   {
      // extract fluid and immersed solid meshes
   }
   feLog("-------------------------------\nFluid Mesh\n");
	m.mesh->PrintCharacteristics();
   feLog("-------------------------------\n");

   // check to ensure at most two materials have been specified!
   if (fem.Materials() > 2)
   {
      feLogError("More than two materials have been specified. At most two materials can be set for the MFEMiFSI solver.");
      return false;
   }

   // get fluid material
   for (int i=0; i<fem.Materials(); ++i)
   {
      FEMaterial* mat = fem.GetMaterial(i);
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
	m.flowsolver = new INS_Solver(m.mesh, (int) m.orderv, (int) m.orderp, m.kinvis, m.fluid_density);

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
   if (m.direct_solver) m.flowsolver->UseDirectSolver();
   m.flowsolver->SetPreconditioner(m.prec_type, m.direct_solver_prec);

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

	// Set the initial conditions.
	for (int i = 0; i < fem.InitialConditions(); ++i)
	{
		MFEMInitialCondition* ic = dynamic_cast<MFEMInitialCondition*>(fem.InitialCondition(i));
		if (ic && ic->IsActive())
		{
			ic->Apply(this);
		}
	}

	// set boundary conditions and loads
	for (int i=0; i<fem.BoundaryConditions(); ++i)
	{
		MFEMBoundaryCondition* bc = dynamic_cast<MFEMBoundaryCondition*>(fem.BoundaryCondition(i));
		if (bc && bc->IsActive())
		{
			FESurface& surf = *bc->GetSurface();
			FEFacetSet* set = surf.GetFacetSet();
			int id = m.mesh->GetSurfaceAttribute(set); assert(id > 0);
			bc->Apply(this, id - 1);
		}
	}
	for (int i = 0; i < fem.ModelLoads(); ++i)
	{
      // nodal loads
      MFEMNodalLoad* nl = dynamic_cast<MFEMNodalLoad*>(fem.ModelLoad(i));
      if (nl && nl->IsActive())
      {
         nl->Apply(this);
      }

      // surface loads
		MFEMLoad* ml = dynamic_cast<MFEMLoad*>(fem.ModelLoad(i));
		if (ml && ml->IsActive())
		{
         FESurface& surf = ml->GetSurface();
         FEFacetSet* set = surf.GetFacetSet();
         int id = m.mesh->GetSurfaceAttribute(set); assert(id > 0);
         ml->Apply(this, id - 1);
		}

      // body loads
      MFEMBodyLoad* bl = dynamic_cast<MFEMBodyLoad*>(fem.ModelLoad(i));
      if (bl && bl->IsActive())
      {
         for (int j=0; j<bl->Domains(); ++j)
         {
            FEDomain* dom = bl->Domain(j);
            int id = m.mesh->GetDomainAttribute(dom); assert(id > 0);
            bl->Apply(this, id - 1);
         }         
      }
	}

   // Zero mean pressure if active
   if (m.zero_mean_pressure) m.flowsolver->AddZeroMeanPressureCondition();

   // Fluid grid functions
	m.u_gf = m.flowsolver->GetCurrentVelocity();
	m.p_gf = m.flowsolver->GetCurrentPressure();

   // Add immersed solid
   if (m.immersed_solid)
   {
      // create solid mesh
      // TODO

      // create solid material
      MFEMSolidMaterial* solid_mat = nullptr;
      for (int i=0; i<fem.Materials(); ++i)
      {
         FEMaterial* fe_mat = fem.GetMaterial(i);
         solid_mat = dynamic_cast<MFEMSolidMaterial*>(fe_mat);
         if (solid_mat)
         {
            break;
         }
      }
      if (!solid_mat)
      {
         feLogError("No immersed solid material found");
         return false;
      }

      // determine which material
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

      // add immersed solid to the flow solver
      FEScalarCoefficient *init_rho_coeff = solid_mat->GetDensityCoeff();
      m.flowsolver->AddImmersedSolid(m.solid_mesh,m.orderd,init_rho_coeff,m.mat_type,m.form_type,m.subtract_viscous_stress);
   
      // immersed solid grid functions
      m.d_gf = m.flowsolver->GetImmersedSolidDisplacement();
      m.vib_gf = m.flowsolver->GetImmersedSolidVelocity();
   }   

	// setup output in paraview format
   std::string fname = "results_mfem_solver";
   if (m.immersed_solid)
   {
      std::string fluid_filename = fname + "/fluid";
      m.pvdc = new ParaViewDataCollection(fluid_filename, m.mesh);
   }
   else
   {
      m.pvdc = new ParaViewDataCollection(fname, m.mesh);
   }	
	m.pvdc->SetDataFormat(VTKFormat::BINARY32);
	m.pvdc->SetHighOrderOutput(true);
	m.pvdc->SetLevelsOfDetail(m.orderv);
	m.pvdc->SetCycle(0);
	m.pvdc->SetTime(0);
	m.pvdc->RegisterField("Velocity", m.u_gf);
	m.pvdc->RegisterField("Pressure", m.p_gf);
   // TODO: add more outputs (curlv, wss, others)
	m.pvdc->Save();

   if (m.immersed_solid)
   {
      std::string solid_filename = fname + "/solid";
      m.sdc = new ParaViewDataCollection(solid_filename, m.solid_mesh);  
      m.sdc->SetDataFormat(VTKFormat::BINARY);
      m.sdc->SetHighOrderOutput(true);
      m.sdc->SetLevelsOfDetail(m.orderd);
      m.sdc->RegisterField("Displacement", m.d_gf);
      m.sdc->RegisterField("Velocity", m.vib_gf);
      m.sdc->SetCycle(0);
      m.sdc->SetTime(0);
      m.sdc->Save();
   }

	double dt = fem.GetTime().timeIncrement;
	m.flowsolver->Setup(dt);
   
	return true;
}

bool MFEMSolver::SolveStep()
{
	// Get the FEBio model and mesh objects
	FEModel& fem = *GetFEModel();
	FEMesh& femesh = fem.GetMesh();

	FEAnalysis* analysis = fem.GetCurrentStep();
	int timeStep = analysis->m_ntimesteps;

	double t = fem.GetTime().currentTime;
	double dt = fem.GetTime().timeIncrement;
	double t0 = t - dt;

	// note that this function needs the previous time point and will update t0.
	// After this call it should be that t0 == t
	m.flowsolver->Step(t0, dt, timeStep);
   
	// save data
	m.pvdc->SetCycle(timeStep + 1);
	m.pvdc->SetTime(t0);
	m.pvdc->Save();
   if (m.immersed_solid)
   {
      m.sdc->SetCycle(timeStep + 1);
      m.sdc->SetTime(t0);
      m.sdc->Save();
   }

	// copy data to FEBio mesh
	Vector& p = *m.p_gf;
	Vector& v = *m.u_gf;
	int dofP = fem.GetDOFIndex("p");
	int dofVX = fem.GetDOFIndex("vx");
	int dofVY = fem.GetDOFIndex("vy");
	int dofVZ = fem.GetDOFIndex("vz");
   int dofUX;
   int dofUY;
   int dofUZ;
   int dofVIMX;
   int dofVIMY;
   int dofVIMZ;
   Vector& d = *m.d_gf;
   Vector& vib = *m.vib_gf;
   if(m.immersed_solid)
   {
      dofUX = fem.GetDOFIndex("im_ux");
      dofUY = fem.GetDOFIndex("im_uy");
      dofUZ = fem.GetDOFIndex("im_uz");
      dofVIMX = fem.GetDOFIndex("im_vx");
      dofVIMY = fem.GetDOFIndex("im_vy");
      dofVIMZ = fem.GetDOFIndex("im_vz");
   }
	int NN = femesh.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		FENode& node = femesh.Node(i);
		node.set(dofP, p[i]);
		node.set(dofVX, v[i]);
		node.set(dofVY, v[i + NN]);
		node.set(dofVZ, v[i + 2*NN]);
      if (m.immersed_solid)
      {
         node.set(dofUX, d[i]);
         node.set(dofUY, d[i + NN]);
         node.set(dofUZ, d[i + 2 * NN]);
         node.set(dofVIMX, v[i]);
         node.set(dofVIMY, v[i + NN]);
         node.set(dofVIMZ, v[i + 2 * NN]);
      }
	}
	fem.Update();

	feLog("%11s %11s\n", "Time", "dt");
	feLog("%.5E %.5E\n", t, dt);

	// print out relative L2 error at final time step
	if (timeStep == analysis->m_ntime - 1)
	{
		m.flowsolver->PrintTimingData();
	}

	return true;
}
