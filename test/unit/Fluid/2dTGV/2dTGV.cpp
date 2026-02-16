// 2D Taylor Green Vortex Problem

#include "MFEMFSI/IncompressibleFluid/IncompressibleFluidSolver.hpp"

using namespace mfem;
using namespace ifs_solver;
using namespace std;

struct s_NavierContext
{
   double kinvis = 0.025;
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
   double dt = 1E-6;
   double t_final = 1E-5; // 10*dt

   // mesh sizes
   // int ser_ref_levels = 2; // 1/2^3
   // int ser_ref_levels = 3; // 1/2^4
   int ser_ref_levels = 4; // 1/2^5   
   // int ser_ref_levels = 5; // 1/2^6

   bool visualize = false;
} ctx;

void vel_exact(const Vector &x, double t, Vector &u)
{
   double x1 = x(0);
   double x2 = x(1);

   u(0) = -1.0*sin(2*M_PI*x2)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   u(1) = sin(2*M_PI*x1)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
}

void grad_vel_exact(const Vector &x, double t, DenseMatrix &gradu)
{
   double x1 = x(0);
   double x2 = x(1);

   gradu(0,0) = 0.0;
   gradu(0,1) = -1.0*2.0*M_PI*cos(2*M_PI*x2)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   gradu(1,0) = 2.0*M_PI*cos(2*M_PI*x1)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   gradu(1,1) = 0.0;
}

void grad_vel_divstress_exact(const Vector &x, double t, DenseMatrix &gradu)
{
   double x1 = x(0);
   double x2 = x(1);

   gradu(0,0) = 0.0;
   gradu(0,1) = -1.0*2.0*M_PI*cos(2*M_PI*x2)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   gradu(1,0) = 2.0*M_PI*cos(2*M_PI*x1)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   gradu(1,1) = 0.0;

   DenseMatrix graduT(gradu);
   graduT.Transpose();

   // gradu + graduT
   gradu += graduT;

   // nu *(gradu + graduT)
   gradu *= ctx.kinvis;   
}

double p_exact(const Vector &x, double t)
{
   double x1 = x(0);
   double x2 = x(1);
   return -cos(2.0*M_PI*x1)*cos(2.0*M_PI*x2)*exp(-8.0*ctx.kinvis*M_PI*M_PI*t);
}

void grad_vel_dot_n_exact(const Vector &x, double t, Vector &gradun)
{
   double x1 = x(0);
   double x2 = x(1);

   DenseMatrix gradu(x.Size());
   gradu(0,0) = 0.0;
   gradu(0,1) = -1.0*2.0*M_PI*cos(2*M_PI*x2)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   gradu(1,0) = 2.0*M_PI*cos(2*M_PI*x1)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   gradu(1,1) = 0.0;

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

   gradu.Mult(nor,gradun);

   // -pn
   Vector pn(x.Size());
   pn = nor;
   pn *= p_exact(x,t);
   gradun.Add(-1.0,pn);
}

void grad_vel_divstress_dot_n_exact(const Vector &x, double t, Vector &gradun)
{
   double x1 = x(0);
   double x2 = x(1);
   gradun.SetSize(x.Size());
   gradun = 0.0;

   DenseMatrix gradu(x.Size());
   gradu(0,0) = 0.0;
   gradu(0,1) = -1.0*2.0*M_PI*cos(2*M_PI*x2)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   gradu(1,0) = 2.0*M_PI*cos(2*M_PI*x1)*exp(-4.0*ctx.kinvis*M_PI*M_PI*t);
   gradu(1,1) = 0.0;

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

   gradu.Mult(nor,gradun);

   // -pn
   Vector pn(x.Size());
   pn = nor;
   pn *= p_exact(x,t);
   gradun.Add(-1.0,pn);
}

void zero_grad(const Vector &x, double t, DenseMatrix &gradu)
{
   double x1 = x(0);
   double x2 = x(1);

   gradu.SetSize(x.Size());
   gradu = 0.0;
}

void zero_vec(const Vector &x, double t, Vector &u)
{
   u.SetSize(x.Size());
   u = 0.0;
}

double ComputeVectorH1Error(GridFunction *u_gf, VectorCoefficient *ex_u, MatrixCoefficient *ex_gradu, const IntegrationRule *irs[]);

void VisualizeField(socketstream &sock, const char *vishost, int visport,
   ParGridFunction &gf, const char *title,
   int x, int y, int w, int h,
   const char *keys = NULL, bool vec = false);

int main(int argc, char *argv[])
{   
   Mpi::Init(argc, argv);
   Hypre::Init();

   Mesh *mesh = new Mesh("tgv.msh");
   for (int i = 0; i < ctx.ser_ref_levels; ++i)
   {
      mesh->UniformRefinement();
   }

   if (Mpi::Root())
   {
      mesh->PrintCharacteristics();
   }

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
   Array<int> attr(pmesh->bdr_attributes.Max());
   attr = 1;
   flowsolver->AddTractionBC(grad_vel_divstress_dot_n_exact,attr);

   double t = 0.0;
   double dt = ctx.dt;
   double t_final = ctx.t_final;
   bool last_step = false;

   // setup restart files
   int start_step = 1;
   flowsolver->Setup(dt);

   // setup output data
   ParGridFunction *u_gf = flowsolver->GetCurrentVelocity();
   ParGridFunction *p_gf = flowsolver->GetCurrentPressure();

   std::string fname;
   if (ctx.equal_order)
   {
      fname = "EO_results";
   }
   else
   {
      fname = "Mixed_results";
   }
   std::stringstream ss;
   ss << std::fixed << std::setprecision(1) << ctx.rhoinf;
   fname += "_GenAlpha_" + ss.str();
   std::setprecision(8);

   if (ctx.prec_type == ifs_solver::PreconditionerType::SIMPLE)
   {
      fname += "_SIMPLE";
   }
   else if (ctx.prec_type == ifs_solver::PreconditionerType::SCR)
   {
      fname += "_SCR";
   }
   else if (ctx.prec_type == ifs_solver::PreconditionerType::aSCR)
   {
      fname += "_aSCR";
   }
   else if (ctx.prec_type == ifs_solver::PreconditionerType::ILU)
   {
      fname += "_ILU";
   }
   else if (ctx.prec_type == ifs_solver::PreconditionerType::Jacobi)
   {
      fname += "_Jacobi";
   }
   else
   {
      fname += "_noPrec";
   }

   // setup output in paraview format
   ParaViewDataCollection pvdc(fname, pmesh);
   pvdc.SetDataFormat(VTKFormat::BINARY32);
   pvdc.SetHighOrderOutput(true);
   pvdc.SetLevelsOfDetail(ctx.order);
   pvdc.SetCycle(0);
   pvdc.SetTime(t);
   pvdc.RegisterField("Velocity", u_gf);
   pvdc.RegisterField("Pressure", p_gf);
   pvdc.Save();

   // setup visualization
   socketstream sout, sout2, sout3, vis_s, vis_us;
   char vishost[] = "localhost";
   int  visport   = 19916;
   int Wx = 0, Wy = 0; // window position
   int Ww = 500, Wh = 500; // window size
   if (ctx.visualize)
   {
      // plot
      VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                  Wx, Wy, Ww, Wh, "Rcamj");  
      MPI_Barrier(pmesh->GetComm());
      VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                  Wx+500, Wy, Ww, Wh, "Rcamj"); 
      MPI_Barrier(pmesh->GetComm()); 
   }

   // time stepping
   for (int step = start_step; !last_step; ++step)
   {
      if (t + dt >= t_final - 0.5 * dt)
      {
         last_step = true;
      }

      t+= dt;
      flowsolver->Step(t,dt,step);

      // save data
      pvdc.SetCycle(step);
      pvdc.SetTime(t);
      pvdc.Save();

      if (ctx.visualize)
      {
         // // plot
         VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                     Wx, Wy, Ww, Wh, "Rcamj");  
         MPI_Barrier(pmesh->GetComm());
         VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                     Wx+500, Wy, Ww, Wh, "Rcamj"); 
         MPI_Barrier(pmesh->GetComm()); 
      }
   }

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

   // compute velocity h1 norm
   MatrixFunctionCoefficient gradu_excoeff(pmesh->Dimension(),grad_vel_exact);
   MatrixFunctionCoefficient zerograd_excoeff(pmesh->Dimension(),zero_grad);
   VectorFunctionCoefficient zerovec_excoeff(pmesh->Dimension(),zero_vec);

   // set time on coefficients
   gradu_excoeff.SetTime(t);
   zerograd_excoeff.SetTime(t);
   zerovec_excoeff.SetTime(t);

   double h1_err_u = ComputeVectorH1Error(u_gf,&u_excoeff,&gradu_excoeff,irs);
   double h1_norm_u = ComputeVectorH1Error(u_gf,&zerovec_excoeff,&zerograd_excoeff,irs);

   if (Mpi::Root())
   {
      std::cout << setprecision(16);
      // relative L2 norm
      std::cout << "|| u_h - u_ex || / || u_ex || = " << err_u / norm_u << "\n";
      std::cout << "|| p_h - p_ex || / || p_ex || = " << err_p / norm_p << "\n";

      // h1 norm for velocity
      std::cout << "|| u_h - u_ex ||_h1 / || u_ex ||_h1 = " << h1_err_u / h1_norm_u << "\n";

      // absolute L2 error
      std::cout << "|| u_h - u_ex || = " << err_u << "\n";
      std::cout << "|| p_h - p_ex || = " << err_p << "\n";
   }

   flowsolver->PrintTimingData();

   // delete pmesh;
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