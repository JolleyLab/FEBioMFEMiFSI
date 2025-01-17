#include "MFEMSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include "ins_solver.hpp"
#include "MFEMMesh.h"
#include "MFEMInitialCondition.h"
#include "MFEMBoundaryCondition.h"
#include "MFEMLoad.h"
#include <FECore/log.h>

using namespace ins_solver_serial;
using namespace std;
using namespace std::placeholders;
using namespace mfem;

class MFEMSolver::Imp
{
public:
	static Imp* This;
	Imp() { This = this; }

public:
	MFEMMesh* mesh = nullptr;
	INS_Solver* flowsolver = nullptr;
	GridFunction* p_gf = nullptr;
	GridFunction* u_gf = nullptr;
	ParaViewDataCollection* pvdc = nullptr;

public:
	int nonlinear_form = INS_Solver::NonlinearFormType::Convective;
	bool div_stress = true;
	
	bool total_integ = false; // This was true for the 2D case
	bool equal_order = true;
	bool supg_flag = true;
	bool pspg_flag = true; // equal order must be true too!
	bool vms_flag = true;
	bool mass_stab = true;

	bool scr_flag = false;
	bool simple_flag = true;
	bool ascr_flag = false;

	// material parameters?
	double fluid_density = 1.0;
	double dyn_viscosity = 10.0;
	double kinvis = dyn_viscosity / fluid_density; // = 0.025 for 2d Case
	double reynolds = 1.0 / kinvis;
	double rhoinf = 0.5;

	// solver tolerances
	double atol = 1e-12;
	double rtol = 1e-8;

	// prec tolerances
	double atol_prec = 1e-12;
	double rtol_prec = 1e-8;

	// equal order
	int order = 1;

	int ts_method = INS_Solver::TS_Method::GeneralizedAlpha;
};

MFEMSolver::Imp* MFEMSolver::Imp::This = nullptr;

// TODO: Add solver parameters
BEGIN_FECORE_CLASS(MFEMSolver, FESolver)
	ADD_PARAMETER(m.atol, "atol");
	ADD_PARAMETER(m.rtol, "rtol");
END_FECORE_CLASS();

MFEMSolver::MFEMSolver(FEModel* fem) : FESolver(fem), m(*(new MFEMSolver::Imp))
{

}

MFEMSolver::~MFEMSolver()
{
	delete m.mesh;
	delete m.flowsolver;
	delete m.pvdc;
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

void MFEMSolver::AddTractionBC(mfem::VectorCoefficient* coeff, mfem::Array<int>& attr)
{
	m.flowsolver->AddTractionBC(coeff, attr);
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
	FEMesh& femesh = fem.GetMesh();
	m.mesh = CreateMFEMMesh3D(femesh);
	if (m.mesh == nullptr) return false;
	m.mesh->PrintCharacteristics();

	// create the mfem solver class
	if (m.equal_order)
	{
		m.flowsolver = new INS_Solver(m.mesh, m.order, m.order, m.kinvis);
	}
	else
	{
		m.flowsolver = new INS_Solver(m.mesh, m.order, m.order - 1, m.kinvis);
	}

	// set time stepping method
	if (m.ts_method == INS_Solver::TS_Method::GeneralizedAlpha)
	{
		m.flowsolver->UseGenAlpha(m.rhoinf);
	}
	else
	{
		feLogError("Just Gen Alpha for serial solver");
		return false;
	}

	// set tolerances
	m.flowsolver->SetLinearSolveAbsTol(m.atol);
	m.flowsolver->SetLinearSolveRelTol(m.rtol);
	m.flowsolver->SetLinearSolveMaxIt(200);

	m.flowsolver->SetNewtonSolveAbsTol(m.atol);
	m.flowsolver->SetNewtonSolveRelTol(m.rtol);
	m.flowsolver->SetNewtonSolveMinIt(3);
	m.flowsolver->SetNewtonSolveMaxIt(20);

	// NOTE: This was not used in 2d example. Add solver parameter for this.
	m.flowsolver->UseSIMPLEPreconditioner();

	// set nonlinear form type
	m.flowsolver->SetNonlinearFormType(m.nonlinear_form);

	if (m.total_integ) m.flowsolver->UseTotalIntegrator();

	// stabilizations
	if (m.supg_flag) m.flowsolver->UseSUPG();
	if (m.pspg_flag)
	{
		if (!m.equal_order)
		{
			feLogError("Error: pspg needs equal order elements!");
			return false;
		}
		m.flowsolver->UsePSPG();
	}
	if (m.vms_flag) m.flowsolver->UseVMS();
	if (m.mass_stab) m.flowsolver->UseMassStabilization();

	// set viscous term type
	if (m.div_stress)
	{
		m.flowsolver->SetViscousTermType(INS_Solver::ViscousTermType::SymGradientForm);
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
		MFEMLoad* ml = dynamic_cast<MFEMLoad*>(fem.ModelLoad(i));
		if (ml && ml->IsActive())
		{
			FESurface& surf = ml->GetSurface();
			FEFacetSet* set = surf.GetFacetSet();
			int id = m.mesh->GetSurfaceAttribute(set); assert(id > 0);
			ml->Apply(this, id - 1);
		}
	}

	m.u_gf = m.flowsolver->GetCurrentVelocity();
	m.p_gf = m.flowsolver->GetCurrentPressure();

	std::string fname;
	if (m.equal_order)
	{
		fname = "EO_results";
	}
	else
	{
		fname = "Mixed_results";
	}
	if (m.ts_method == INS_Solver::TS_Method::GeneralizedAlpha)
	{
		std::stringstream ss;
		ss << std::fixed << std::setprecision(1) << m.rhoinf;
		fname += "_GenAlpha_" + ss.str();
		std::setprecision(8);
	}
	else
	{
		fname += "_BDF2";
	}

	if (m.simple_flag)
	{
		fname += "_SIMPLE";
	}
	else if (m.scr_flag)
	{
		fname += "_SCR";
	}
	else if (m.ascr_flag)
	{
		fname += "_aSCR";
	}
	else
	{
		MFEM_ABORT("Error: unknown preconditioner specificed");
	}

	// setup output in paraview format
	m.pvdc = new ParaViewDataCollection(fname, m.mesh);
	m.pvdc->SetDataFormat(VTKFormat::BINARY32);
	m.pvdc->SetHighOrderOutput(true);
	m.pvdc->SetLevelsOfDetail(m.order);
	m.pvdc->SetCycle(0);
	m.pvdc->SetTime(0);
	m.pvdc->RegisterField("Velocity", m.u_gf);
	m.pvdc->RegisterField("Pressure", m.p_gf);
	m.pvdc->Save();

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

	// copy data to FEBio mesh
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
		node.set(dofVX, v[i]);
		node.set(dofVY, v[i + NN]);
		node.set(dofVZ, v[i + 2*NN]);
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
