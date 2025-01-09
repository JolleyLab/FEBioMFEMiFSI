#include "MFEMSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FEDomain.h>
#include <FECore/FEAnalysis.h>
#include "ins_solver.hpp"
#include "MFEMMesh.h"
#include <math.h>
#include <FECore/log.h>
#include <functional>

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
	Mesh* mesh = nullptr;
	INS_Solver* flowsolver = nullptr;
	GridFunction* p_gf = nullptr;
	GridFunction* u_gf = nullptr;
	VectorFunctionCoefficient* u_excoeff = nullptr;
	FunctionCoefficient* p_excoeff = nullptr;
	ParaViewDataCollection* pvdc = nullptr;

public:
	double kinvis = 0.025;
	double reynolds = 1.0 / kinvis;
	double rhoinf = 0.5;

	int nonlinear_form = INS_Solver::NonlinearFormType::Convective;
	bool div_stress = true;
	
	bool total_integ = true;
	bool equal_order = true;
	bool supg_flag = equal_order;
	bool pspg_flag = equal_order; // equal order must be true too!
	bool vms_flag = equal_order;
	bool mass_stab = equal_order;

	bool scr_flag = false;
	bool simple_flag = true;
	bool ascr_flag = false;

	// solver tolerances
	double atol = 1e-10;
	double rtol = 1e-6;

	// prec tolerances
	double atol_prec = 1e-10;
	double rtol_prec = 1e-6;

	// equal order
	int order = 1;

	int ts_method = INS_Solver::TS_Method::GeneralizedAlpha;

public:
	void vel_exact(const Vector& x, double t, Vector& u)
	{
		double x1 = x(0);
		double x2 = x(1);

		u(0) = -1.0 * sin(2 * M_PI * x2) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		u(1) = sin(2 * M_PI * x1) * exp(-4.0 * kinvis * M_PI * M_PI * t);
	}

	void grad_vel_exact(const Vector& x, double t, DenseMatrix& gradu)
	{
		double x1 = x(0);
		double x2 = x(1);

		gradu(0, 0) = 0.0;
		gradu(0, 1) = -1.0 * 2.0 * M_PI * cos(2 * M_PI * x2) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		gradu(1, 0) = 2.0 * M_PI * cos(2 * M_PI * x1) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		gradu(1, 1) = 0.0;
	}

	void grad_vel_divstress_exact(const Vector& x, double t, DenseMatrix& gradu)
	{
		double x1 = x(0);
		double x2 = x(1);

		gradu(0, 0) = 0.0;
		gradu(0, 1) = -1.0 * 2.0 * M_PI * cos(2 * M_PI * x2) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		gradu(1, 0) = 2.0 * M_PI * cos(2 * M_PI * x1) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		gradu(1, 1) = 0.0;

		DenseMatrix graduT(gradu);
		graduT.Transpose();

		// gradu + graduT
		gradu += graduT;

		// nu *(gradu + graduT)
		gradu *= kinvis;
	}

	double p_exact(const Vector& x, double t)
	{
		double x1 = x(0);
		double x2 = x(1);
		return -cos(2.0 * M_PI * x1) * cos(2.0 * M_PI * x2) * exp(-8.0 * kinvis * M_PI * M_PI * t);
	}

	void grad_vel_dot_n_exact(const Vector& x, double t, Vector& gradun)
	{
		double x1 = x(0);
		double x2 = x(1);

		DenseMatrix gradu(x.Size());
		gradu(0, 0) = 0.0;
		gradu(0, 1) = -1.0 * 2.0 * M_PI * cos(2 * M_PI * x2) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		gradu(1, 0) = 2.0 * M_PI * cos(2 * M_PI * x1) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		gradu(1, 1) = 0.0;

		gradu *= kinvis;

		double eps = 1e-8;
		Vector nor(x.Size());
		nor = 0.0;
		if (abs(x2 - (-0.5)) < eps)
		{
			nor(1) = -1.0;
		}
		else if (abs(x2 - (0.5)) < eps)
		{
			nor(1) = 1.0;
		}
		else if (abs(x1 - (0.5)) < eps)
		{
			nor(0) = 1.0;
		}
		else if (abs(x1 - (-0.5)) < eps)
		{
			nor(0) = -1.0;
		}

		gradu.Mult(nor, gradun);

		// -pn
		Vector pn(x.Size());
		pn = nor;
		pn *= p_exact(x, t);
		gradun.Add(-1.0, pn);
	}

	// NOTE: The AddTractionBC doesn't support functionals, so we can't bind this function. 
	//       Instead, we'll use a different mechanism.
	static void grad_vel_divstress_dot_n_exact(const Vector& x, double t, Vector& gradun)
	{
		This->Grad_vel_divstress_dot_n_exact(x, t, gradun);
	}

	void Grad_vel_divstress_dot_n_exact(const Vector& x, double t, Vector& gradun)
	{
		double x1 = x(0);
		double x2 = x(1);
		gradun.SetSize(x.Size());
		gradun = 0.0;

		DenseMatrix gradu(x.Size());
		gradu(0, 0) = 0.0;
		gradu(0, 1) = -1.0 * 2.0 * M_PI * cos(2 * M_PI * x2) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		gradu(1, 0) = 2.0 * M_PI * cos(2 * M_PI * x1) * exp(-4.0 * kinvis * M_PI * M_PI * t);
		gradu(1, 1) = 0.0;

		DenseMatrix graduT(gradu);
		graduT.Transpose();

		// gradu + graduT
		gradu += graduT;

		// nu *(gradu + graduT)
		gradu *= kinvis;

		double eps = 1e-8;
		Vector nor(x.Size());
		nor = 0.0;
		if (abs(x2 - (-0.5)) < eps)
		{
			nor(1) = -1.0;
		}
		else if (abs(x2 - (0.5)) < eps)
		{
			nor(1) = 1.0;
		}
		else if (abs(x1 - (0.5)) < eps)
		{
			nor(0) = 1.0;
		}
		else if (abs(x1 - (-0.5)) < eps)
		{
			nor(0) = -1.0;
		}

		gradu.Mult(nor, gradun);

		// -pn
		Vector pn(x.Size());
		pn = nor;
		pn *= p_exact(x, t);
		gradun.Add(-1.0, pn);
	}
};

MFEMSolver::Imp* MFEMSolver::Imp::This = nullptr;

// TODO: Add solver parameters
BEGIN_FECORE_CLASS(MFEMSolver, FESolver)
	ADD_PARAMETER(m.atol, "atol");
	ADD_PARAMETER(m.rtol, "rtol");
END_FECORE_CLASS();

void zero_grad(const Vector& x, double t, DenseMatrix& gradu)
{
	double x1 = x(0);
	double x2 = x(1);

	gradu.SetSize(x.Size());
	gradu = 0.0;
}

void zero_vec(const Vector& x, double t, Vector& u)
{
	u.SetSize(x.Size());
	u = 0.0;
}

double ComputeVectorH1Error(GridFunction* u_gf, VectorCoefficient* ex_u, MatrixCoefficient* ex_gradu, const IntegrationRule* irs[]);

void ComputeRelativeMass_Momentum_Error(GridFunction& v_gf, double fluid_density, VectorCoefficient& vel_exact_coeff, Vector& result);

double ComputeVectorH1Error(GridFunction* u_gf, VectorCoefficient* ex_u, MatrixCoefficient* ex_gradu, const IntegrationRule* irs[])
{
	double error_temp = 0.0;
	const FiniteElement* fe;
	ElementTransformation* Tr;
	DenseMatrix grad;
	int intorder;
	int dim = u_gf->FESpace()->GetMesh()->SpaceDimension();
	DenseMatrix mat(dim);
	Vector uvec(dim), vec(dim);

	for (int i = 0; i < u_gf->FESpace()->GetNE(); i++)
	{
		fe = u_gf->FESpace()->GetFE(i);
		Tr = u_gf->FESpace()->GetElementTransformation(i);
		intorder = 2 * fe->GetOrder() + 3; // <--------
		const IntegrationRule* ir;
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
			const IntegrationPoint& ip = ir->IntPoint(j);
			Tr->SetIntPoint(&ip);
			u_gf->GetVectorGradient(*Tr, grad);
			ex_gradu->Eval(mat, *Tr, ip);
			mat.Add(-1.0, grad);

			u_gf->GetVectorValue(*Tr, ip, vec);
			ex_u->Eval(uvec, *Tr, ip);
			vec -= uvec;

			error_temp += ip.weight * Tr->Weight() * (mat.FNorm2() + vec * vec);
		}
	}

	return sqrt(abs(error_temp));
}

void ComputeRelativeMass_Momentum_Error(GridFunction& v_gf, double fluid_density, VectorCoefficient& vel_exact_coeff, Vector& result)
{
	const FiniteElement* fe;
	ElementTransformation* Tr;
	DenseMatrix gradv;
	int intorder;
	int dim = v_gf.FESpace()->GetMesh()->SpaceDimension();
	double div_v;
	double div_int = 0.0;
	double mass_int = 0.0;
	double momentum_exact_int = 0.0;
	double momentum_int = 0.0;
	double momentum_diff = 0.0;
	double div_err = 0.0;

	double div_err_l2 = 0.0;

	Vector v(dim), v_exact(dim);

	gradv.SetSize(dim);
	gradv = 0.0;

	double p = 1.0;

	for (int i = 0; i < v_gf.FESpace()->GetNE(); i++)
	{
		fe = v_gf.FESpace()->GetFE(i);
		Tr = v_gf.FESpace()->GetElementTransformation(i);
		intorder = 2 * fe->GetOrder() + 3; // <--------
		const IntegrationRule* ir;
		ir = &(IntRules.Get(fe->GetGeomType(), intorder));
		for (int j = 0; j < ir->GetNPoints(); j++)
		{
			const IntegrationPoint& ip = ir->IntPoint(j);
			Tr->SetIntPoint(&ip);
			v_gf.GetVectorGradient(*Tr, gradv);

			vel_exact_coeff.Eval(v_exact, *Tr, ip);
			v_gf.GetVectorValue(*Tr, ip, v);

			// compute div v
			div_v = gradv.Trace();

			// div u
			div_int += ip.weight * Tr->Weight() * div_v;

			// compute mass
			mass_int += ip.weight * Tr->Weight() * fluid_density;

			// compute momentum_exact
			momentum_exact_int += ip.weight * Tr->Weight() * fluid_density * v_exact.Norml2();

			// momentum error
			v -= v_exact;
			double loc_error = v.Norml2();
			momentum_diff += ip.weight * Tr->Weight() * loc_error;
		}
	}

	// compute relative errors
	result.SetSize(2);
	// glob mass in this case = 1
	result(0) = div_int / 1.0;
	result(1) = momentum_diff / momentum_exact_int;
}

MFEMSolver::MFEMSolver(FEModel* fem) : FESolver(fem), m(*(new MFEMSolver::Imp))
{

}

MFEMSolver::~MFEMSolver()
{
	delete m.mesh;
	delete m.flowsolver;
	delete m.u_excoeff;
	delete m.p_excoeff;
	delete m.pvdc;
}

bool MFEMSolver::Init()
{
	// create the mfem mesh from the FEBio mesh
	FEModel& fem = *GetFEModel();
	FEMesh& femesh = fem.GetMesh();
	m.mesh = CreateMFEMMesh(femesh);
	if (m.mesh == nullptr) return false;

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

	// Set the initial condition.
	GridFunction* u_ic = m.flowsolver->GetCurrentVelocity();
	m.u_excoeff = new VectorFunctionCoefficient(m.mesh->Dimension(), bind(&Imp::vel_exact, &m, _1, _2, _3));
	m.u_excoeff->SetTime(0.0);
	u_ic->ProjectCoefficient(*m.u_excoeff);

	GridFunction* p_ic = m.flowsolver->GetCurrentPressure();
	m.p_excoeff = new FunctionCoefficient((bind(&Imp::p_exact, &m, _1, _2)));
	m.p_excoeff->SetTime(0.0);
	p_ic->ProjectCoefficient(*m.p_excoeff);

	// set boundary conditions
	Array<int> attr(m.mesh->bdr_attributes.Max());
	attr = 1;
	m.flowsolver->AddTractionBC(&Imp::grad_vel_divstress_dot_n_exact, attr);

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
		node.set(dofVZ, 0);
	}
	fem.Update();

	printf("%11s %11s\n", "Time", "dt");
	printf("%.5E %.5E\n", t, dt);
	fflush(stdout);

	// print out relative L2 error at final time step
	if (timeStep == analysis->m_ntime - 1)
	{
		double err_u = 0.0;
		double err_p = 0.0;
		m.u_excoeff->SetTime(t);
		m.p_excoeff->SetTime(t);

		// compute relative L2 error on finer mesh
		int order_err = m.order + 3;
		const IntegrationRule* irs[Geometry::NumGeom];
		for (int i = 0; i < Geometry::NumGeom; ++i)
		{
			irs[i] = &(IntRules.Get(i, order_err));
		}
		err_u = m.u_gf->ComputeL2Error(*m.u_excoeff, irs);
		double norm_u = ComputeLpNorm(2., *m.u_excoeff, *m.mesh, irs);
		err_p = m.p_gf->ComputeL2Error(*m.p_excoeff, irs);
		double norm_p = ComputeLpNorm(2., *m.p_excoeff, *m.mesh, irs);

		// compute velocity h1 norm
		MatrixFunctionCoefficient gradu_excoeff(m.mesh->Dimension(), bind(&Imp::grad_vel_exact, &m, _1, _2, _3));
		MatrixFunctionCoefficient zerograd_excoeff(m.mesh->Dimension(), zero_grad);
		VectorFunctionCoefficient zerovec_excoeff(m.mesh->Dimension(), zero_vec);

		// set time on coefficients
		gradu_excoeff.SetTime(t);
		zerograd_excoeff.SetTime(t);
		zerovec_excoeff.SetTime(t);

		double h1_err_u = ComputeVectorH1Error(m.u_gf, &*m.u_excoeff, &gradu_excoeff, irs);
		double h1_norm_u = ComputeVectorH1Error(m.u_gf, &zerovec_excoeff, &zerograd_excoeff, irs);

		std::cout << setprecision(16);
		// relative L2 norm
		std::cout << "|| u_h - u_ex || / || u_ex || = " << err_u / norm_u << "\n";
		std::cout << "|| p_h - p_ex || / || p_ex || = " << err_p / norm_p << "\n";

		// h1 norm for velocity
		std::cout << "|| u_h - u_ex ||_h1 / || u_ex ||_h1 = " << h1_err_u / h1_norm_u << "\n";

		// absolute L2 error
		std::cout << "|| u_h - u_ex || = " << err_u << "\n";
		std::cout << "|| p_h - p_ex || = " << err_p << "\n";

		m.flowsolver->PrintTimingData();
	}

	return true;
}
