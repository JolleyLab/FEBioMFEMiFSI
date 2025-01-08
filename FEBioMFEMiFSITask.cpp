#include "FEBioMFEMiFSITask.h"

// Navier 2D Taylor Green Vortex Problem
#include "ins_solver.hpp"
#include <fstream>
#include <math.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FEEdge.h>

using namespace mfem;
using namespace ins_solver_serial;
using namespace std;

struct s_NavierContext
{
	// original
	double kinvis = 0.025;
	double reynolds = 1.0 / kinvis;

	int ts_method = INS_Solver::TS_Method::GeneralizedAlpha;
	// // double rhoinf = 0.0; // use for space converge
	double rhoinf = 0.5;
	// // double rhoinf = 1.0; // use for time converge

	int nonlinear_form = INS_Solver::NonlinearFormType::Convective;
	bool div_stress = true;
	bool visualize = false;
	bool tet_mesh = false;

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


	// spatial convergence //
	// implicit method for all spatial convergence!
	double dt = 1E-6;
	double t_final = 1E-5; // 10*dt

	// int ser_ref_levels = 0; // 1/2^1

	// int ser_ref_levels = 1; // 1/2^2

	// for stabilizations do the 4 below

	// int ser_ref_levels = 2; // 1/2^3
	// int ser_ref_levels = 3; // 1/2^4
	int ser_ref_levels = 4; // 1/2^5   
	// int ser_ref_levels = 5; // 1/2^6


	// int ser_ref_levels = 6; // 1/2^7

	// equal order
	int order = 1;
	// int order = 2;
	// int order = 3;

	// mixed
	// int order = 2;
	// int order = 3;
	// int order = 4;



	// // // temporal convergence //
	// int order = 5;
	// int ser_ref_levels = 4; // 1/2^5
	// bool equal_order = false;

	// // fixed dt
	// // logspace from 10^-2 to 10^-3 with 4 points
	// double t_final = 2E-2;
	// // double dt = 0.01;
	// // double dt = 0.005;
	// // double dt = 0.0025;
	// double dt = 0.00125;   
} ctx;

void vel_exact(const Vector& x, double t, Vector& u)
{
	double x1 = x(0);
	double x2 = x(1);

	u(0) = -1.0 * sin(2 * M_PI * x2) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	u(1) = sin(2 * M_PI * x1) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
}

void grad_vel_exact(const Vector& x, double t, DenseMatrix& gradu)
{
	double x1 = x(0);
	double x2 = x(1);

	gradu(0, 0) = 0.0;
	gradu(0, 1) = -1.0 * 2.0 * M_PI * cos(2 * M_PI * x2) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	gradu(1, 0) = 2.0 * M_PI * cos(2 * M_PI * x1) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	gradu(1, 1) = 0.0;
}

void grad_vel_divstress_exact(const Vector& x, double t, DenseMatrix& gradu)
{
	double x1 = x(0);
	double x2 = x(1);

	gradu(0, 0) = 0.0;
	gradu(0, 1) = -1.0 * 2.0 * M_PI * cos(2 * M_PI * x2) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	gradu(1, 0) = 2.0 * M_PI * cos(2 * M_PI * x1) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	gradu(1, 1) = 0.0;

	DenseMatrix graduT(gradu);
	graduT.Transpose();

	// gradu + graduT
	gradu += graduT;

	// nu *(gradu + graduT)
	gradu *= ctx.kinvis;
}

double p_exact(const Vector& x, double t)
{
	double x1 = x(0);
	double x2 = x(1);
	return -cos(2.0 * M_PI * x1) * cos(2.0 * M_PI * x2) * exp(-8.0 * ctx.kinvis * M_PI * M_PI * t);
}

void grad_vel_dot_n_exact(const Vector& x, double t, Vector& gradun)
{
	double x1 = x(0);
	double x2 = x(1);

	DenseMatrix gradu(x.Size());
	gradu(0, 0) = 0.0;
	gradu(0, 1) = -1.0 * 2.0 * M_PI * cos(2 * M_PI * x2) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	gradu(1, 0) = 2.0 * M_PI * cos(2 * M_PI * x1) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	gradu(1, 1) = 0.0;

	gradu *= ctx.kinvis;

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

void grad_vel_divstress_dot_n_exact(const Vector& x, double t, Vector& gradun)
{
	double x1 = x(0);
	double x2 = x(1);
	gradun.SetSize(x.Size());
	gradun = 0.0;

	DenseMatrix gradu(x.Size());
	gradu(0, 0) = 0.0;
	gradu(0, 1) = -1.0 * 2.0 * M_PI * cos(2 * M_PI * x2) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	gradu(1, 0) = 2.0 * M_PI * cos(2 * M_PI * x1) * exp(-4.0 * ctx.kinvis * M_PI * M_PI * t);
	gradu(1, 1) = 0.0;

	DenseMatrix graduT(gradu);
	graduT.Transpose();

	// gradu + graduT
	gradu += graduT;

	// nu *(gradu + graduT)
	gradu *= ctx.kinvis;


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

bool FEBioMFEMiFSITask::Init(const char* szfile)
{
	return true;
}

// This function constructs a mfem mesh from the FEBioMesh.
// Note that this currently assumes a 2D geometry! 
// TODO: Generalize this for 3D meshes.
Mesh* CreateMesh(FEMesh& febioMesh)
{
	int nodes = febioMesh.Nodes();
	int elems = febioMesh.Elements();

	// In FEBio, it is not necessary to provide the entire boundary. 
	// I don't know if this will be an issue for mfem. 
	int bdrElems = 0;
	for (int i = 0; i < febioMesh.SegmentSets(); ++i)
	{
		bdrElems += febioMesh.SegmentSet(i).Segments();
	}

	Mesh* mesh = new Mesh(2, nodes, elems, bdrElems, 2);
	for (int i = 0; i < nodes; ++i)
	{
		FENode& node = febioMesh.Node(i);
		vec3d& pos = node.m_r0;
		double coords[3] = { pos.x, pos.y, pos.z };
		mesh->AddVertex(coords);
	}

	int indices[FEElement::MAX_NODES] = { 0 };
	for (int i = 0; i < febioMesh.SegmentSets(); ++i)
	{
		FESegmentSet& set = febioMesh.SegmentSet(i);
		for (int j = 0; j < set.Segments(); ++j)
		{
			FESegmentSet::SEGMENT& line = set.Segment(j);
			indices[0] = line.node[0];
			indices[1] = line.node[1];
			mesh->AddBdrSegment(indices[0], indices[1], i + 1);
		}
	}

	for (int n = 0; n < febioMesh.Domains(); ++n)
	{
		FEDomain& dom = febioMesh.Domain(n);
		int elemsInDomain = dom.Elements();
		for (int i = 0; i < elemsInDomain; i++) {
			FEElement& el = dom.ElementRef(i);
			for (int j = 0; j < el.Nodes(); j++) {
				indices[j] = el.m_node[j];
			}

			switch (el.Shape())
			{
			case ET_QUAD4: mesh->AddQuad(indices[0], indices[1], indices[2], indices[3], n + 1); break;
			default:
				assert(false);
			}
		}
	}

	mesh->FinalizeTopology(false);
	mesh->Finalize(1, true);

	return mesh;
}

//! Run the task.
bool FEBioMFEMiFSITask::Run()
{
	// Get the FEBio model and mesh objects
	FEModel& fem = *GetFEModel();
	FEMesh& femesh = fem.GetMesh();

	// construct an mfem mesh from the febio mesh
	Mesh* mesh = CreateMesh(femesh);

	// do some refinement
	for (int i = 0; i < ctx.ser_ref_levels; ++i)
	{
		mesh->UniformRefinement();
	}

	// Create the flow solver.
	INS_Solver* flowsolver = nullptr;
	if (ctx.equal_order)
	{
		flowsolver = new INS_Solver(mesh, ctx.order, ctx.order, ctx.kinvis);
	}
	else
	{
		flowsolver = new INS_Solver(mesh, ctx.order, ctx.order - 1, ctx.kinvis);
	}

	// set time stepping method
	if (ctx.ts_method == INS_Solver::TS_Method::GeneralizedAlpha)
	{
		flowsolver->UseGenAlpha(ctx.rhoinf);
	}
	else
	{
		MFEM_ABORT("Just Gen Alpha for serial solver");
	}

	// set tolerances
	flowsolver->SetLinearSolveAbsTol(ctx.atol);
	flowsolver->SetLinearSolveRelTol(ctx.rtol);
	flowsolver->SetLinearSolveMaxIt(200);

	flowsolver->SetNewtonSolveAbsTol(ctx.atol);
	flowsolver->SetNewtonSolveRelTol(ctx.rtol);
	flowsolver->SetNewtonSolveMinIt(3);
	flowsolver->SetNewtonSolveMaxIt(20);

	// // preconditioner tolerances
	// flowsolver->SetPrecAINVAbsTol(ctx.atol_prec);
	// flowsolver->SetPrecAINVRelTol(ctx.rtol_prec);
	// flowsolver->SetPrecSINVAbsTol(ctx.atol_prec);
	// flowsolver->SetPrecSINVRelTol(ctx.rtol_prec);

	// // Set block preconditioner
	// if (ctx.scr_flag) flowsolver->UseSCRPreconditioner();
	// else if (ctx.simple_flag) flowsolver->UseSIMPLEPreconditioner();  
	// else if (ctx.ascr_flag) flowsolver->Use_aSCRPreconditioner();
	// else MFEM_ABORT("Error");

	// set nonlinear form type
	flowsolver->SetNonlinearFormType(ctx.nonlinear_form);

	if (ctx.total_integ) flowsolver->UseTotalIntegrator();

	// stabilizations
	if (ctx.supg_flag) flowsolver->UseSUPG();
	if (ctx.pspg_flag)
	{
		MFEM_ASSERT(ctx.equal_order, "Error: pspg needs equal order elements!");
		flowsolver->UsePSPG();
	}
	if (ctx.vms_flag) flowsolver->UseVMS();
	if (ctx.mass_stab) flowsolver->UseMassStabilization();

	// set viscous term type
	if (ctx.div_stress)
	{
		flowsolver->SetViscousTermType(INS_Solver::ViscousTermType::SymGradientForm);
	}

	// Set the initial condition.
	GridFunction* u_ic = flowsolver->GetCurrentVelocity();
	VectorFunctionCoefficient u_excoeff(mesh->Dimension(), vel_exact);
	u_excoeff.SetTime(0.0);
	u_ic->ProjectCoefficient(u_excoeff);
	GridFunction* p_ic = flowsolver->GetCurrentPressure();
	FunctionCoefficient p_excoeff(p_exact);
	p_excoeff.SetTime(0.0);
	p_ic->ProjectCoefficient(p_excoeff);

	// set boundary conditions
	Array<int> attr(mesh->bdr_attributes.Max());
	attr = 1;
	flowsolver->AddTractionBC(grad_vel_divstress_dot_n_exact, attr);

	double t = 0.0;
	double dt = ctx.dt;
	double t_final = ctx.t_final;
	bool last_step = false;

	flowsolver->Setup(dt);

	GridFunction* u_gf = flowsolver->GetCurrentVelocity();
	GridFunction* p_gf = flowsolver->GetCurrentPressure();

	std::string fname;
	if (ctx.equal_order)
	{
		fname = "EO_results";
	}
	else
	{
		fname = "Mixed_results";
	}
	if (ctx.ts_method == INS_Solver::TS_Method::GeneralizedAlpha)
	{
		std::stringstream ss;
		ss << std::fixed << std::setprecision(1) << ctx.rhoinf;
		fname += "_GenAlpha_" + ss.str();
		std::setprecision(8);
	}
	else
	{
		fname += "_BDF2";
	}

	if (ctx.simple_flag)
	{
		fname += "_SIMPLE";
	}
	else if (ctx.scr_flag)
	{
		fname += "_SCR";
	}
	else if (ctx.ascr_flag)
	{
		fname += "_aSCR";
	}
	else
	{
		MFEM_ABORT("Error: unknown preconditioner specificed");
	}

	// setup output in paraview format
	ParaViewDataCollection pvdc(fname, mesh);
	pvdc.SetDataFormat(VTKFormat::BINARY32);
	pvdc.SetHighOrderOutput(true);
	pvdc.SetLevelsOfDetail(ctx.order);
	pvdc.SetCycle(0);
	pvdc.SetTime(t);
	pvdc.RegisterField("Velocity", u_gf);
	pvdc.RegisterField("Pressure", p_gf);
	pvdc.Save();

	// time stepping
	for (int step = 0; !last_step; ++step)
	{
		if (t + dt >= t_final - 0.5 * dt)
		{
			last_step = true;
		}

		flowsolver->Step(t, dt, step);

		// save data
		pvdc.SetCycle(step + 1);
		pvdc.SetTime(t);
		pvdc.Save();

		printf("%11s %11s\n", "Time", "dt");
		printf("%.5E %.5E\n", t, dt);
		fflush(stdout);
	}

	// print out relative L2 error at final time step
	double err_u = 0.0;
	double err_p = 0.0;
	u_excoeff.SetTime(t);
	p_excoeff.SetTime(t);

	// compute relative L2 error on finer mesh
	int order_err = ctx.order + 3;
	const IntegrationRule* irs[Geometry::NumGeom];
	for (int i = 0; i < Geometry::NumGeom; ++i)
	{
		irs[i] = &(IntRules.Get(i, order_err));
	}
	err_u = u_gf->ComputeL2Error(u_excoeff, irs);
	double norm_u = ComputeLpNorm(2., u_excoeff, *mesh, irs);
	err_p = p_gf->ComputeL2Error(p_excoeff, irs);
	double norm_p = ComputeLpNorm(2., p_excoeff, *mesh, irs);

	// compute velocity h1 norm
	MatrixFunctionCoefficient gradu_excoeff(mesh->Dimension(), grad_vel_exact);
	MatrixFunctionCoefficient zerograd_excoeff(mesh->Dimension(), zero_grad);
	VectorFunctionCoefficient zerovec_excoeff(mesh->Dimension(), zero_vec);

	// set time on coefficients
	gradu_excoeff.SetTime(t);
	zerograd_excoeff.SetTime(t);
	zerovec_excoeff.SetTime(t);

	double h1_err_u = ComputeVectorH1Error(u_gf, &u_excoeff, &gradu_excoeff, irs);
	double h1_norm_u = ComputeVectorH1Error(u_gf, &zerovec_excoeff, &zerograd_excoeff, irs);

	std::cout << setprecision(16);
	// relative L2 norm
	std::cout << "|| u_h - u_ex || / || u_ex || = " << err_u / norm_u << "\n";
	std::cout << "|| p_h - p_ex || / || p_ex || = " << err_p / norm_p << "\n";

	// h1 norm for velocity
	std::cout << "|| u_h - u_ex ||_h1 / || u_ex ||_h1 = " << h1_err_u / h1_norm_u << "\n";

	// absolute L2 error
	std::cout << "|| u_h - u_ex || = " << err_u << "\n";
	std::cout << "|| p_h - p_ex || = " << err_p << "\n";

	flowsolver->PrintTimingData();

	// delete pmesh;
	delete mesh;
	delete flowsolver;
	return true;
}
