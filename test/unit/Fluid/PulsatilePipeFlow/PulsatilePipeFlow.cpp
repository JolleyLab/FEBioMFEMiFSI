// 3D pipe flow with windkessel BC at outflow

#include "MFEMFSI/IncompressibleFluid/IncompressibleFluidSolver.hpp"

using namespace mfem;
using namespace ifs_solver;
using namespace std;

struct IFS_Settings
{
   // blood properties
   double fluid_density = 1.06; // g/cm^2
   double dyn_viscosity = 0.04; // g/cm/s
   double kinvis = dyn_viscosity / fluid_density;

   // pipe geometry
   double D = 4.0; // cm
   double Radius = D/2.0; //cm
   double L = 30.0; // 30cm pipe

   // flow rate and windkessel parameters
   double Q0 = 10.0;
   double Rp = 0.1;
   double Rd = 1.0;
   double C = 1.0/4.0/M_PI;
   double p0 = 0.0;
   double pd = 0.0;
   double vmax = Q0*M_PI*2.0/Radius/Radius;
   double cfl_max = 1.0;

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
   std::string device_config = "cpu";

   // solver tolerances
   double atol = 1e-6;
   double rtol = 1e-3;
   int minIt = 2;
   int maxIt = 10;

   int prec_type = ifs_solver::PreconditionerType::SIMPLE;
   double lin_atol = 1e-6;
   double lin_rtol = 1e-3;
   int lin_maxIt = 200;

   double prec_Ainv_atol = 1e-6;
   double prec_Ainv_rtol = 1e-3;
   double prec_Sinv_atol = 1e-6;
   double prec_Sinv_rtol = 1e-3;
   
   // time stepping
   int ts_method = ifs_solver::TimeSteppingMethod::GeneralizedAlpha;
   double rhoinf = 0.5;
   bool set_by_cfl = true;
   double t_final = 3.0;
   double dt = 5E-3; // dt is not used, its set by CFL = 1.0
   int step_save = 10;

   // mesh sizes
   enum class MeshResolution { Coarse, Medium, Fine };
   MeshResolution mesh_resolution = MeshResolution::Coarse;

   bool visualization = false;
} ctx;

double flow_rate(double t)
{
   return ctx.Q0 * sin(t/ctx.Rd/ctx.C/2.0) * sin(t/ctx.Rd/ctx.C/2.0);
}

void parabolic_velo(const Vector &x, double t, Vector &u)
{
   double xi = x(0);
   double yi = x(1);
   double zi = x(2);

   u.SetSize(x.Size());
   u = 0.0;

   // parabolic profile (vmax = 1) on circle of radius R has flow rate = pi*r*r/2
   double flow_rate_normalization = M_PI*ctx.Radius*ctx.Radius/2.0;

   // parabolic profile
   double r = sqrt(xi*xi + yi*yi); // pipe along z direction
   u(2) = (flow_rate(t) / flow_rate_normalization) * (1.0 - ((r*r) / (ctx.Radius * ctx.Radius))); // z-velocity
}

void VisualizeField(socketstream &sock, const char *vishost, int visport,
   ParGridFunction &gf, const char *title,
   int x, int y, int w, int h,
   const char *keys = NULL, bool vec = false)
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

int main(int argc, char *argv[])
{   
   Mpi::Init(argc, argv);
   Hypre::Init();

   // Device configuration
   Device device(ctx.device_config.c_str());
   if (Mpi::Root()) 
   {
      mfem::out << "Using device: " << ctx.device_config << std::endl;
      device.Print();
      mfem::out << std::endl;
   }

   std::string mesh_file;
   switch (ctx.mesh_resolution)
   {
      case IFS_Settings::MeshResolution::Coarse:
         mesh_file = "pipe_coarse.msh";
         break;
      case IFS_Settings::MeshResolution::Medium:
         mesh_file = "pipe_medium.msh";
         break;
      case IFS_Settings::MeshResolution::Fine:
         mesh_file = "pipe_fine.msh";
         break;
      default:
         MFEM_ABORT("Unknown mesh resolution");
   }
   
   int inlet = 1;
   int outlet = 2;
   int wall = 3;

   Mesh *mesh = new Mesh(mesh_file);
   if (Mpi::Root())
   {
      mesh->PrintCharacteristics();
   }

   double hmin,hmax,kappa_min,kappa_max;
   mesh->GetCharacteristics(hmin,hmax,kappa_min,kappa_max);

   // compute dt based on CFL num
   // no need to worry about time constant for windkessel (its > 0.05)
   if (ctx.set_by_cfl) ctx.dt = hmin * ctx.cfl_max / ctx.vmax;

   // partition the mesh
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);

   // Create the flow solver  
   IncompressibleFluidSolver *flowsolver = new IncompressibleFluidSolver(pmesh, ctx.order, ctx.equal_order ? ctx.order : ctx.order - 1, ctx.kinvis, ctx.fluid_density, device, ctx.orderingtype);

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

   // Set the initial condition. (0 everywhere)
   
   // init pressure func for windkessel
   GridFunction p_init(*flowsolver->GetCurrentPressure()); // copy constructor
   GridFunctionCoefficient p_init_func(&p_init);

   // set boundary conditions
   // Add inflow condition
   Array<int> inflow_attr(pmesh->bdr_attributes.Max());
   inflow_attr = 0;
   inflow_attr[inlet-1] = 1; // attr 1
   flowsolver->AddVelDirichletBC(parabolic_velo,inflow_attr);
   
   // add bc at outlet
   Array<int> attr_outflow(pmesh->bdr_attributes.Max());
   attr_outflow = 0;
   attr_outflow[outlet-1] = 1; // attr 2

   // windkessel bc
   // init flow rate for windkessel
   double init_flow_rate = flowsolver->ComputeMassFlowThruSurf(*flowsolver->GetCurrentVelocity(),attr_outflow,false,ctx.fluid_density,false);
   flowsolver->AddWindkesselBC(attr_outflow,ctx.Rp,ctx.Rd,ctx.C,ctx.p0,ctx.pd,init_flow_rate);
   flowsolver->UseBackflowStabilization(attr_outflow);
   
   // add no slip bc
   Array<int> no_slip_attr(pmesh->bdr_attributes.Max());
   no_slip_attr = 0;
   no_slip_attr[wall-1] = 1; // attr 3
   flowsolver->AddNoSlipBC(no_slip_attr);

   // delete serial mesh
   delete mesh;

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
   flowsolver->ComputeWSS(*u_gf,*p_gf);
   ParGridFunction *wss_gf = flowsolver->GetCurrentWSS();

   string fname = std::string("results_pipe_flow_np") + std::to_string(Mpi::WorldSize());
   ParaViewDataCollectionBase *dc = nullptr;
   #ifdef MFEM_USE_HDF5
   dc = new ParaViewHDFDataCollection(fname, pmesh);
   #else
   dc = new ParaViewDataCollection(fname, pmesh);
   #endif
   dc->SetDataFormat(VTKFormat::BINARY32);
   dc->SetLevelsOfDetail(ctx.order);
   dc->SetCycle(0);
   dc->SetTime(t);
   dc->RegisterField("Velocity", u_gf);
   dc->RegisterField("Pressure", p_gf);
   dc->RegisterField("WSS", wss_gf);
   dc->Save();

   // setup visualization
   socketstream sout, sout2, sout3;
   char vishost[] = "localhost";
   int  visport   = 19916;
   int Wx = 0, Wy = 0; // window position
   int Ww = 500, Wh = 500; // window size

   if (ctx.visualization)
   {
      // plot
      VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                  Wx, Wy, Ww, Wh, "campppppppppppppppppppppppp",false);
      MPI_Barrier(pmesh->GetComm()); 
      VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                  Wx+500, Wy, Ww, Wh, "campppppppppppppppppppppppp",false);
      MPI_Barrier(pmesh->GetComm()); 
      VisualizeField(sout3, vishost, visport, *wss_gf, "WSS",
         Wx+1000, Wy, Ww, Wh,"campppppppppppppppppppppppp"); 
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
      flowsolver->Step(t,dt,step,false);

      if ((step) % ctx.step_save == 0 || last_step)
      {
         flowsolver->ComputeWSS(*u_gf,*p_gf);
         dc->SetCycle(step);
         dc->SetTime(t);
         dc->Save();
      }

      if (ctx.visualization)
      {
         // plot
         VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                     Wx, Wy, Ww, Wh, "campppppppppppppppppppppppp",false);
         MPI_Barrier(pmesh->GetComm()); 
         VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                     Wx+500, Wy, Ww, Wh, "campppppppppppppppppppppppp",false);
         MPI_Barrier(pmesh->GetComm()); 
         flowsolver->ComputeWSS(*u_gf,*p_gf);
         VisualizeField(sout3, vishost, visport, *wss_gf, "WSS",
            Wx+1000, Wy, Ww, Wh,"campppppppppppppppppppppppp"); 
         MPI_Barrier(pmesh->GetComm());
      }
   }

   flowsolver->PrintTimingData();

   delete dc;
   delete pmesh;
   delete flowsolver;

   return 0;
}