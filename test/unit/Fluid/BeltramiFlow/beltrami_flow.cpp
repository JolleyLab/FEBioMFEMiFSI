// 3D Beltrami flow

#include "MFEMFSI/IncompressibleFluid/IncompressibleFluidSolver.hpp"

using namespace mfem;
using namespace ifs_solver;
using namespace std;

struct s_NavierContext
{
   double fluid_density = 1.0;
   double dyn_viscosity = 10.0;
   double kinvis = dyn_viscosity / fluid_density;
   double reynolds = 1.0 / kinvis;

   // VMS
   int stab_type = ifs_solver::StabilizationType::VMS;
   int order = 1;
   bool equal_order = true;

   // // Stable 
   // int stab_type = ifs_solver::StabilizationType::NONE;
   // int order = 2;
   // bool equal_order = false;

   int nonlinear_form = ifs_solver::NonlinearFormType::Convective;
   int viscous_term_type = ifs_solver::ViscousTermType::StressForm; // {LaplaceForm = 0, StressForm = 1}
   int orderingtype = Ordering::byVDIM;

   // solver tolerances
   double atol = 1e-10;
   double rtol = 1e-6;
   int minIt = 1;
   int maxIt = 10;

   int prec_type = ifs_solver::PreconditionerType::SIMPLE;
   double lin_atol = 1e-10;
   double lin_rtol = 1e-6;
   int lin_maxIt = 200;

   double prec_Ainv_atol = 1e-6;
   double prec_Ainv_rtol = 1e-3;
   double prec_Sinv_atol = 1e-6;
   double prec_Sinv_rtol = 1e-3;
   
   // time stepping
   int ts_method = ifs_solver::TimeSteppingMethod::GeneralizedAlpha;
   double rhoinf = 0.5;
   double dt = 10.0E-6;
   double t_final = 10.0E-5; // 10*dt

   // mesh sizes
   // int nel_1d = pow(2,2); // h = 2*pi / 2^2
   // int nel_1d = pow(2,3); // h = 2*pi / 2^3
   int nel_1d = pow(2,4); // h = 2*pi / 2^4
   // int nel_1d = pow(2,5); // h = 2*pi / 2^5
   // int nel_1d = pow(2,6); // h = 2*pi / 2^6
} ctx;

void vel_exact(const Vector &xval, double t, Vector &u)
{
   double x = xval(0);
   double y = xval(1);
   double z = xval(2);

   double c = M_PI*exp(-M_PI*M_PI*ctx.kinvis*t/4.0);

   u(0) = -1.0*c/4.0 * ( exp(M_PI*x/4.0)*sin(M_PI*(y + 2.0*z)/4.0) + exp(M_PI*z/4.0)*cos(M_PI*(x + 2.0*y)/4.0) );
   u(1) = -1.0*c/4.0 * ( exp(M_PI*y/4.0)*sin(M_PI*(2.0*x + z)/4.0) + exp(M_PI*x/4.0)*cos(M_PI*(y + 2.0*z)/4.0) );
   u(2) = -1.0*c/4.0 * ( exp(M_PI*z/4.0)*sin(M_PI*(x + 2.0*y)/4.0) + exp(M_PI*y/4.0)*cos(M_PI*(2.0*x + z)/4.0) );
}

double p_exact(const Vector &xval, double t)
{
   double x = xval(0);
   double y = xval(1);
   double z = xval(2);

   double time_part = -ctx.fluid_density*M_PI*M_PI/32.0*exp(-M_PI*M_PI*t*ctx.kinvis/2.0);

   double space_part = exp(M_PI*x/2.0) + exp(M_PI*y/2.0) + exp(M_PI*z/2.0)
                        + 2.0*exp(M_PI*(x + y)/4.0)*cos(M_PI*(y + 2.0*z)/4.0)*sin(M_PI*(2.0*x + z)/4.0)
                        + 2.0*exp(M_PI*(x + z)/4.0)*cos(M_PI*(x + 2.0*y)/4.0)*sin(M_PI*(y + 2.0*z)/4.0)
                        + 2.0*exp(M_PI*(y + z)/4.0)*cos(M_PI*(2.0*x + z)/4.0)*sin(M_PI*(x + 2.0*y)/4.0);


   return time_part*space_part;
}

void traction_bc(const Vector &xval, double t, Vector &tract)
{
   double x = xval(0);
   double y = xval(1);
   double z = xval(2);

   tract.SetSize(xval.Size());
   tract = 0.0;

   double p = p_exact(xval,t);

   tract(0) = p + (M_PI*M_PI*ctx.dyn_viscosity*exp(-(M_PI*M_PI*ctx.kinvis*t)/4.0)*(exp((M_PI*x)/4.0)*sin((M_PI*(y + 2.0*z))/4.0) - exp((M_PI*z)/4.0)*sin((M_PI*(x + 2.0*y))/4.0)))/8.0;
   tract(1) = 2.0*ctx.dyn_viscosity*((M_PI*exp(-(M_PI*M_PI*ctx.kinvis*t)/4.0)*((M_PI*exp((M_PI*x)/4.0)*cos((M_PI*(y + 2.0*z))/4.0))/4.0 + (M_PI*exp((M_PI*y)/4.0)*cos((M_PI*(2.0*x + z))/4.0))/2.0))/8.0 - (M_PI*exp(-(M_PI*M_PI*ctx.kinvis*t)/4.0)*((M_PI*exp((M_PI*z)/4.0)*sin((M_PI*(x + 2.0*y))/4.0))/2.0 - (M_PI*exp((M_PI*x)/4.0)*cos((M_PI*(y + 2.0*z))/4.0))/4.0))/8.0);
   tract(2) = 2.0*ctx.dyn_viscosity*((M_PI*exp(-(M_PI*M_PI*ctx.kinvis*t)/4.0)*((M_PI*exp((M_PI*x)/4.0)*cos((M_PI*(y + 2.0*z))/4.0))/2.0 + (M_PI*exp((M_PI*z)/4.0)*cos((M_PI*(x + 2.0*y))/4.0))/4.0))/8.0 - (M_PI*exp(-(M_PI*M_PI*ctx.kinvis*t)/4.0)*((M_PI*exp((M_PI*y)/4.0)*sin((M_PI*(2.0*x + z))/4.0))/2.0 - (M_PI*exp((M_PI*z)/4.0)*cos((M_PI*(x + 2.0*y))/4.0))/4.0))/8.0);
}

void VisualizeField(socketstream &sock, const char *vishost, int visport,
   ParGridFunction &gf, const char *title,
   int x, int y, int w, int h,
   const char *keys = NULL, bool vec = false);

int main(int argc, char *argv[])
{   
   Mpi::Init(argc, argv);
   Hypre::Init();

   Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(ctx.nel_1d, ctx.nel_1d, ctx.nel_1d, Element::HEXAHEDRON, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI));
   // Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(ctx.nel_1d, ctx.nel_1d, ctx.nel_1d, Element::TETRAHEDRON, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI));
   // Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(ctx.nel_1d, ctx.nel_1d, ctx.nel_1d, Element::WEDGE, 2.0*M_PI, 2.0*M_PI, 2.0*M_PI));
   if (Mpi::Root()) mesh->PrintCharacteristics();

   // partition the mesh
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh; // delete serial mesh

   // Create the flow solver.
   IncompressibleFluidSolver *flowsolver = nullptr;    
   if (ctx.equal_order)
   {
      flowsolver = new IncompressibleFluidSolver(pmesh, ctx.order, ctx.order, ctx.kinvis, ctx.orderingtype);
   }
   else
   {
      flowsolver = new IncompressibleFluidSolver(pmesh, ctx.order, ctx.order-1, ctx.kinvis, ctx.orderingtype);
   }

   // Time stepping method
   flowsolver->SetTimeSteppingMethod(ctx.ts_method);
	flowsolver->SetRhoInf(ctx.rhoinf);

	// set tolerances
	flowsolver->SetLinearSolveAbsTol(ctx.lin_atol);
	flowsolver->SetLinearSolveRelTol(ctx.lin_rtol);
	flowsolver->SetLinearSolveMaxIt(ctx.lin_maxIt);

	flowsolver->SetNewtonSolveAbsTol(ctx.atol);
	flowsolver->SetNewtonSolveRelTol(ctx.rtol);
	flowsolver->SetNewtonSolveMinIt(ctx.minIt);
	flowsolver->SetNewtonSolveMaxIt(ctx.maxIt);

   flowsolver->SetPrecAINVAbsTol(ctx.prec_Ainv_atol);
   flowsolver->SetPrecAINVRelTol(ctx.prec_Ainv_rtol);
   flowsolver->SetPrecSINVAbsTol(ctx.prec_Sinv_atol);
   flowsolver->SetPrecSINVRelTol(ctx.prec_Sinv_rtol);

	// set linear solver and preconditioner
   flowsolver->SetPreconditioner(ctx.prec_type, false);

	// set nonlinear form type
	flowsolver->SetNonlinearFormType(ctx.nonlinear_form);

   // use total integrator
	flowsolver->UseTotalIntegrator();

	// stabilizations
   switch (ctx.stab_type) 
   {
      case ifs_solver::StabilizationType::NONE:
         break;

      case ifs_solver::StabilizationType::VMS:
      default:
         flowsolver->UseVMS();
         break;
   }

	// set viscous term type, default is StressForm
   flowsolver->SetViscousTermType(ctx.viscous_term_type);

   // Set the initial condition.
   ParGridFunction *u_ic = flowsolver->GetCurrentVelocity();
   VectorFunctionCoefficient u_excoeff(pmesh->Dimension(), vel_exact);
   u_excoeff.SetTime(0.0);
   u_ic->ProjectCoefficient(u_excoeff);
   ParGridFunction *p_ic = flowsolver->GetCurrentPressure();
   FunctionCoefficient p_excoeff(p_exact);
   p_excoeff.SetTime(0.0);
   p_ic->ProjectCoefficient(p_excoeff);

   // set boundary conditions
   // Add dirichlet velocity bc
   Array<int> vel_bc_attr(pmesh->bdr_attributes.Max());
   vel_bc_attr = 1;
   vel_bc_attr[4] = 0;
   flowsolver->AddVelDirichletBC(vel_exact, vel_bc_attr);

   // neumann bc
   Array<int> vel_n_attr(pmesh->bdr_attributes.Max());
   vel_n_attr = 0;
   vel_n_attr[4] = 1;

   flowsolver->AddTractionBC(traction_bc,vel_n_attr);

   double t = 0.0;
   double dt = ctx.dt;
   double t_final = ctx.t_final;
   bool last_step = false;

   // setup
   flowsolver->Setup(dt);

   ParGridFunction *u_gf = flowsolver->GetCurrentVelocity();
   ParGridFunction *p_gf = flowsolver->GetCurrentPressure();

   ParaViewDataCollection dc("results_beltrami_flow", pmesh);
   dc.SetDataFormat(VTKFormat::BINARY32);
   dc.SetHighOrderOutput(true);
   dc.SetLevelsOfDetail(ctx.order);
   dc.SetCycle(0);
   dc.SetTime(t);
   dc.RegisterField("velocity", u_gf);
   dc.RegisterField("pressure", p_gf);
   dc.Save();

   MPI_Barrier(pmesh->GetComm()); 

   for (int step = 1; !last_step; ++step)
   {
      if (t + dt >= t_final - 0.5 * dt)
      {
         last_step = true;
      }

      t += dt;
      flowsolver->Step(t,dt,step,false);

      if ((step) % 5 == 0 || last_step)
      {
         dc.SetCycle(step);
         dc.SetTime(t);
         dc.Save();
      }
   }
   MPI_Barrier(pmesh->GetComm()); 

   // print out relative L2 error at final time step
   double err_u = 0.0;
   double err_p = 0.0;
   u_excoeff.SetTime(t);
   p_excoeff.SetTime(t);
   
   // compute relative L2 error on finer mesh
   int order_err = ctx.order + 3;
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i=0; i < Geometry::NumGeom; ++i)
   {
      irs[i] = &(IntRules.Get(i, order_err));
   }
   err_u  = u_gf->ComputeL2Error(u_excoeff, irs);
   double norm_u = ComputeGlobalLpNorm(2., u_excoeff, *pmesh, irs);
   err_p  = p_gf->ComputeL2Error(p_excoeff, irs);
   double norm_p = ComputeGlobalLpNorm(2., p_excoeff, *pmesh, irs);

   if (Mpi::Root())
   {
      std::cout << setprecision(16);
      std::cout << "|| u_h - u_ex || / || u_ex || = " << err_u / norm_u << "\n";
      std::cout << "|| p_h - p_ex || / || p_ex || = " << err_p / norm_p << "\n";
      std::cout << "|| u_h - u_ex || = " << err_u << "\n";
      std::cout << "|| p_h - p_ex || = " << err_p << "\n";
   }

   flowsolver->PrintTimingData();

   delete pmesh;
   delete flowsolver;

   return 0;
}

double ComputeVectorH1Error(GridFunction *u_gf, VectorCoefficient *ex_u, MatrixCoefficient *ex_gradu, const IntegrationRule *irs[])
{
   double error_temp = 0.0;
   const FiniteElement *fe;
   ElementTransformation *Tr;
   DenseMatrix grad;
   int intorder;
   int dim = u_gf->FESpace()->GetMesh()->SpaceDimension();
   DenseMatrix mat(dim);
   Vector uvec(dim), vec(dim);

   for (int i = 0; i < u_gf->FESpace()->GetNE(); i++)
   {
      fe = u_gf->FESpace()->GetFE(i);
      Tr = u_gf->FESpace()->GetElementTransformation(i);
      intorder = 2*fe->GetOrder() + 3; // <--------
      const IntegrationRule *ir;
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
         const IntegrationPoint &ip = ir->IntPoint(j);
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

   // parallelization
   double glob_error = 0.0;
   double local_error = abs(error_temp);
   MPI_Allreduce(&local_error, &glob_error, 1, MPI_DOUBLE, MPI_SUM, dynamic_cast<ParGridFunction*>(u_gf)->ParFESpace()->GetComm());

   return sqrt(glob_error);
}

void VisualizeField(socketstream &sock, const char *vishost, int visport,
   ParGridFunction &gf, const char *title,
   int x, int y, int w, int h, const char *keys, bool vec)
{
   ParMesh &pmesh = *gf.ParFESpace()->GetParMesh();
   MPI_Comm comm = pmesh.GetComm();

   int num_procs, myid;
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   bool newly_opened = false;
   int connection_failed;

   do
   {
   if (myid == 0)
   {
   if (!sock.is_open() || !sock)
   {
   sock.open(vishost, visport);
   sock.precision(8);
   newly_opened = true;
   }
   sock << "solution\n";
   }

   pmesh.PrintAsOne(sock);
   gf.SaveAsOne(sock);

   if (myid == 0 && newly_opened)
   {
   sock << "window_title '" << title << "'\n"
   << "window_geometry "
   << x << " " << y << " " << w << " " << h << "\n";
   if ( keys ) { sock << "keys " << keys << "\n"; }
   else { sock << "keys maaAc"; }
   if ( vec ) { sock << "vvv"; }
   sock << endl;
   }

   if (myid == 0)
   {
   connection_failed = !sock && !newly_opened;
   }
   MPI_Bcast(&connection_failed, 1, MPI_INT, 0, comm);
   }
   while (connection_failed);
}

void ComputeRelativeMass_Momentum_Error(GridFunction &v_gf,double fluid_density, VectorCoefficient &vel_exact_coeff, Vector &result)
{
   const FiniteElement *fe;
   ElementTransformation *Tr;
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
      intorder = 2*fe->GetOrder() + 3; // <--------
      const IntegrationRule *ir;
      ir = &(IntRules.Get(fe->GetGeomType(), intorder));
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir->IntPoint(j);
         Tr->SetIntPoint(&ip);
         v_gf.GetVectorGradient(*Tr,gradv);

         vel_exact_coeff.Eval(v_exact,*Tr,ip);
         v_gf.GetVectorValue(*Tr,ip,v);

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

   // parallelization
   double glob_mass = 0.0;
   double glob_div_v = 0.0;
   MPI_Allreduce(&div_int, &glob_div_v, 1, MPI_DOUBLE, MPI_SUM, dynamic_cast<ParGridFunction*>(&v_gf)->ParFESpace()->GetComm());
   MPI_Allreduce(&mass_int, &glob_mass, 1, MPI_DOUBLE, MPI_SUM, dynamic_cast<ParGridFunction*>(&v_gf)->ParFESpace()->GetComm());

   double glob_m_exact = 0.0;
   double glob_m = 0.0;
   MPI_Allreduce(&momentum_exact_int, &glob_m_exact, 1, MPI_DOUBLE, MPI_SUM, dynamic_cast<ParGridFunction*>(&v_gf)->ParFESpace()->GetComm());
   MPI_Allreduce(&momentum_diff, &glob_m, 1, MPI_DOUBLE, MPI_SUM, dynamic_cast<ParGridFunction*>(&v_gf)->ParFESpace()->GetComm());

   // compute relative errors
   result.SetSize(2);
   result(0) = glob_div_v / 1.0;
   result(1) = glob_m / glob_m_exact;
}