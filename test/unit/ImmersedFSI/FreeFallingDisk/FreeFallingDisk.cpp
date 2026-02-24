// Free falling disk

#include "MFEMFSI/IncompressibleFluid/IncompressibleFluidSolver.hpp"
#include "MFEMFSI/ImmersedFSI/ImmersedSolid.hpp"
#include "MFEMFSI/ImmersedFSI/Materials/materials.hpp"

using namespace mfem;
using namespace ifs_solver;
using namespace std;

struct s_NavierContext
{
   double fluid_density = 1000.0; // kg/m^3
   double dyn_viscosity = 0.1;
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
   double rtol = 1e-3;
   int minIt = 2;
   int maxIt = 10;

   int prec_type = ifs_solver::PreconditionerType::SIMPLE;
   double lin_atol = 1e-10;
   double lin_rtol = 1e-3;
   int lin_maxIt = 200;

   double prec_Ainv_atol = 1e-6;
   double prec_Ainv_rtol = 1e-3;
   double prec_Sinv_atol = 1e-6;
   double prec_Sinv_rtol = 1e-3;

   // time stepping
   int ts_method = ifs_solver::TimeSteppingMethod::GeneralizedAlpha;
   double rhoinf = 0.5;
   
   // fluid mesh
   int dim = 2;
   // coarse
   int problem = 1;
   int nx = 48;
   int ny = 36;

   // // medium
   // int problem = 2;
   // int nx = 96;
   // int ny = 72;

   // // fine
   // int problem = 3;
   // int nx = 192;
   // int ny = 144;

   // // very fine
   // int problem = 4;
   // int nx = 384;
   // int ny = 288;

   // int dim = 3; // one option
} ctx;

struct s_ImmersedSolidContext
{
   int order = ctx.order;

   // solid material properties
   double nu = 0.4;
   double E = 1076.0; // Pa
   double lambda = nu*E / ((1.0 + nu)*(1.0-2.0*nu));
   double mu = E / (2.0*(1.0+nu));
   double solid_density = 1200.0; // kg/m^3

   // penalty parameter
   double kappa = mu/10.0; // penalty value (default value)

   bool almost_rigid = false; // (default value)
   bool subtract_visc_stress = true;

   // int form_type = ImmersedSolid::form_type::TotalLagrangian;
   int form_type = ImmersedSolid::form_type::UpdatedLagrangian;
   int vol_type = NearlyIncompressibleMaterial::volumetric_energy::ST91;
} stx;

struct s_ImmersedFSIContext
{
   double gravity = 9.81; // m/s^2

   double dt = 1e-3;
   int step_save = 10;

   int min_it = 2;
   int max_it = 10;

   // double t_final = 1.5; // s
   double t_final = 0.2; // enough for disk to move a significant distance

   bool visualization = false;
   int orderingtype = Ordering::byVDIM;
   // int orderingtype = Ordering::byNODES;
} itx;

void gravity_func(const Vector &x, double t, Vector &u)
{
   u.SetSize(x.Size());
   u = 0.0;
   u(0) = itx.gravity;
}

double vel_zero(const Vector &x, double t)
{
   return 0.0;
}

void temp_func(const Vector &x, double t, Vector &u)
{
   u.SetSize(x.Size());
   u = 1.0;
}

void VisualizeMesh(socketstream &sock, const char *vishost, int visport,
   ParMesh &pmesh, const char *title,
   int x = 0, int y = 0, int w = 400, int h = 400,
   const char *keys = NULL);
void VisualizeField(socketstream &sock, const char *vishost, int visport,
   ParGridFunction &gf, const char *title,
   int x, int y, int w, int h,
   const char *keys = NULL, bool vec = false);

int main(int argc, char *argv[])
{   
   Mpi::Init(argc, argv);
   Hypre::Init();

   //////////// Initialize and parameters ////////////////

   // Read the mesh from the given mesh file.
   Mesh *fluid_mesh = nullptr;
   Mesh *solid_mesh = nullptr;
   if (ctx.dim == 2)
   {
      fluid_mesh = new Mesh(Mesh::MakeCartesian2D(ctx.nx, ctx.ny, Element::QUADRILATERAL, false, 0.16, 0.04));
      solid_mesh = new Mesh("solid_circle_1820quads.vtk");
   }
   else
   {
      fluid_mesh = new Mesh("fluid_mesh_febio.mesh");
      solid_mesh = new Mesh("solid_mesh_febio.mesh");
   }
   int dim = fluid_mesh->Dimension();

   // ensure nodes for each mesh, this is required for GSLIB
   fluid_mesh->SetCurvature(ctx.order);
   solid_mesh->SetCurvature(stx.order);
   
   if (Mpi::Root())
   {
      std::cout << "fluid mesh characteristics =" << std::endl;
      fluid_mesh->PrintCharacteristics();
      std::cout << "solid mesh characteristics =" << std::endl;
      solid_mesh->PrintCharacteristics();
   }   

   // partition the fluid mesh
   ParMesh *fluid_pmesh = new ParMesh(MPI_COMM_WORLD, *fluid_mesh);
   delete fluid_mesh; // delete serial mesh

   // partition the solid mesh
   ParMesh *solid_pmesh = new ParMesh(MPI_COMM_WORLD, *solid_mesh);
   delete solid_mesh; // delete serial mesh

   // Create the flow solver  
   IncompressibleFluidSolver *flowsolver = new IncompressibleFluidSolver(fluid_pmesh, ctx.order, ctx.equal_order ? ctx.order : ctx.order - 1, ctx.kinvis, ctx.fluid_density, ctx.orderingtype);

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
   
   if (ctx.dim == 3)
   {
      // add bcs for fluid
      // left end = traction free

      // top and bottom walls = slip (no velo in vertical direction)
      Array<int> topandbot_attr(flowsolver->GetMesh()->bdr_attributes.Max());
      topandbot_attr = 0;
      topandbot_attr[1] = 1;
      flowsolver->AddVelDirichletBC(vel_zero,topandbot_attr,2); // uy = 0, ux free

      // right end = no slip
      Array<int> right_attr(flowsolver->GetMesh()->bdr_attributes.Max());
      right_attr = 0;
      right_attr[0] = 1;
      flowsolver->AddNoSlipBC(right_attr);

      // front and back walls = slip (no velo in z direction)
      Array<int> z_attr(flowsolver->GetMesh()->bdr_attributes.Max());
      z_attr = 0;
      z_attr[2] = 1;
      flowsolver->AddVelDirichletBC(vel_zero,z_attr,3); // uz = 0, ux,uy free
   }
   else
   {
      // add bcs for fluid
      // left end = traction free
      Array<int> left_attr(flowsolver->GetMesh()->bdr_attributes.Max());
      left_attr = 0;
      left_attr[3] = 1;
      // flowsolver->AddTractionFreeBC(left_attr);   

      // top and bottom walls = slip (no velo in vertical direction)
      Array<int> topandbot_attr(flowsolver->GetMesh()->bdr_attributes.Max());
      topandbot_attr = 0;
      topandbot_attr[0] = 1;
      topandbot_attr[2] = 1;
      flowsolver->AddVelDirichletBC(vel_zero,topandbot_attr,2); // uy = 0, ux free

      // right end = no slip
      Array<int> right_attr(flowsolver->GetMesh()->bdr_attributes.Max());
      right_attr = 0;
      right_attr[1] = 1;
      flowsolver->AddNoSlipBC(right_attr);
   }   

   // body force of gravity over all fluid (and solid!)
   Array<int> domain_attr(flowsolver->GetMesh()->attributes.Max());
   domain_attr = 1;
   flowsolver->AddBodyForce(gravity_func,domain_attr);

   // solid material
   ConstantCoefficient init_rhos(stx.solid_density);
   // material coeffs
   ConstantCoefficient mu_coeff(stx.mu);
   ConstantCoefficient kappa_coeff(stx.kappa);
   ConstantCoefficient lambda_coeff(stx.lambda);
   
   // add immersed solid
   MaterialType *mat_type = new NearlyIncompressibleNeoHookean(stx.form_type,solid_pmesh->Dimension(),stx.vol_type,kappa_coeff,mu_coeff);
   flowsolver->AddImmersedSolid(solid_pmesh,stx.order,&init_rhos,mat_type,stx.form_type,stx.subtract_visc_stress);

   ///////////////////////////////////////////////////////////////////////////////////
   // setup solvers
   flowsolver->Setup(itx.dt);

   ///////////////////////////////////////////////////////////////////////////////////
   // folder names
   std::string folder_name;
   if (ctx.problem == 1)
   {
      folder_name = "result_freefall_coarse";
   }
   else if (ctx.problem == 2)
   {
      folder_name = "result_freefall_medium";
   }
   else if (ctx.problem == 3)
   {
      folder_name = "result_freefall_fine";
   }   
   else
   {
      folder_name = "result_freefall_veryfine";
   }

   if (ctx.dim == 3) folder_name += "_3D";

   folder_name += "_np" + std::to_string(Mpi::WorldSize());

   ///////////////////////////////////////////////////////////////////////////////////
   // Setup visualization
   ParGridFunction *u_gf = flowsolver->GetCurrentVelocity();
   ParGridFunction *p_gf = flowsolver->GetCurrentPressure();

   std::string fluid_filename;
   fluid_filename = folder_name + "/fluid";
   ParaViewDataCollection *fdc = nullptr;
   GridFunction *f_fsi_fluid = nullptr;
   fdc = new ParaViewDataCollection(fluid_filename,fluid_pmesh);
   fdc->SetDataFormat(VTKFormat::BINARY);
   fdc->SetHighOrderOutput(true);
   fdc->SetLevelsOfDetail(ctx.order);
   fdc->SetPrecision(8);
   fdc->RegisterField("Velocity", u_gf);
   fdc->RegisterField("Pressure", p_gf);
   ParGridFunction w_gf(*u_gf);
   if (ctx.dim == 2) flowsolver->ComputeCurl2D(*u_gf, w_gf);
   else flowsolver->ComputeCurl3D(*u_gf, w_gf);
   fdc->RegisterField("Vorticity", &w_gf);

   fdc->SetCycle(0);
   fdc->SetTime(0.0);
   fdc->Save();

   // visualization
   std::string solid_filename;
   solid_filename = folder_name + "/solid";
   ParGridFunction *d_gf = flowsolver->GetImmersedSolidDisplacement();
   ParGridFunction *vib_gf = flowsolver->GetImmersedSolidVelocity();
   ParGridFunction *a_gf = flowsolver->GetImmersedSolidAcceleration();
   // copy solid mesh for outputting on deformed state
   ParMesh current_solid_pmesh(*solid_pmesh);
   GridFunction temp_gf(d_gf->FESpace());
   temp_gf = 0.0;
   ParaViewDataCollection *sdc = new ParaViewDataCollection(solid_filename,&current_solid_pmesh);  
   sdc->SetDataFormat(VTKFormat::BINARY);
   sdc->SetHighOrderOutput(true);
   sdc->SetLevelsOfDetail(stx.order);
   sdc->SetPrecision(8);
   sdc->RegisterField("Displacement", d_gf);
   sdc->RegisterField("Velocity", vib_gf);
   sdc->RegisterField("Acceleration", a_gf);
   sdc->SetCycle(0);
   sdc->SetTime(0.0);
   sdc->Save();

   ///////////////////////////////////////////////////////////////////////////////////
   socketstream sout, sout2, sout3, sout4;
   char vishost[] = "localhost";
   int  visport   = 19916;
   int Ww = 500, Wh = 500; // window size
   int Wx = 0, Wy = 0; // window position

   if(itx.visualization)
   {
      if (ctx.dim == 2)
      {
         VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                  Wx, Wy, Ww, Wh, "Rcj");           
         VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                     Wx+500, Wy, Ww, Wh, "Rcj");
         VisualizeMesh(sout3, vishost, visport, current_solid_pmesh, "Solid Mesh",
            Wx+1000, Wy, Ww, Wh, "pppppppma");
      }
      else
      {
         VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                  Wx, Wy, Ww, Wh, "Rcppppppppppppppppppppppppppppppppp");           
         VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                     Wx+500, Wy, Ww, Wh, "Rcppppppppppppppppppppppppppppppppp");
         VisualizeField(sout3, vishost, visport, *d_gf, "Displacement",
            Wx+1000, Wy, Ww, Wh, "Rcpppppppppppppppppppppppppppppppppm");
      }
      MPI_Barrier(fluid_pmesh->GetComm()); 
   }

   double t = 0.0;
   double dt = itx.dt;
   double t_final = itx.t_final;
   bool last_step = false;

   for (int step = 1; !last_step; ++step)
   {
      if (t + dt >= t_final - dt / 2)
      {
         last_step = true;
      }

      t += dt;
      flowsolver->Step(t,dt,step,false);

      // save fluid and solid results
      if ((step) % itx.step_save == 0 || last_step)
      {
         // fluid
         if (ctx.dim == 2) flowsolver->ComputeCurl2D(*u_gf, w_gf);
         else flowsolver->ComputeCurl3D(*u_gf, w_gf);

         fdc->SetCycle(step);
         fdc->SetTime(t);
         fdc->Save();

         if(sdc)
         {
            // deform current solid mesh
            auto nodes = solid_pmesh->GetNodes();
            if(stx.form_type == ImmersedSolid::form_type::TotalLagrangian)
            {
               add(1.0,*nodes,1.0,*d_gf,temp_gf);
               current_solid_pmesh.SetNodes(temp_gf);
            }
            else
            {
               current_solid_pmesh.SetNodes(*nodes);
            }

            sdc->SetCycle(step);
            sdc->SetTime(t);
            sdc->Save();
         }         
      }

      if(itx.visualization && (step) % itx.step_save == 0)
      {
         if (ctx.dim == 2)
         {
            VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                     Wx, Wy, Ww, Wh, "Rcj");           
            VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                        Wx+500, Wy, Ww, Wh, "Rcj");
            VisualizeMesh(sout3, vishost, visport, current_solid_pmesh, "Solid Mesh",
               Wx+1000, Wy, Ww, Wh, "pppppppma");
         }
         else
         {
            VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                     Wx, Wy, Ww, Wh, "Rcppppppppppppppppppppppppppppppppp");           
            VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                        Wx+500, Wy, Ww, Wh, "Rcppppppppppppppppppppppppppppppppp");
            VisualizeField(sout3, vishost, visport, *d_gf, "Displacement",
                           Wx+1000, Wy, Ww, Wh, "Rcpppppppppppppppppppppppppppppppppm");
         }
         MPI_Barrier(fluid_pmesh->GetComm());
      }      
   }

   flowsolver->PrintTimingData();

   // free memory
   delete flowsolver;
   delete fluid_pmesh;

   delete mat_type;
   delete solid_pmesh;

   delete fdc;
   delete sdc; 

   return 0;
}

void VisualizeMesh(socketstream &sock, const char *vishost, int visport,
                   Mesh &mesh, const char *title,
                   int x, int y, int w, int h, const char * keys)
{
   bool newly_opened = false;
   int connection_failed;

   do
   {
      if (!sock.is_open() || !sock)
      {
         sock.open(vishost, visport);
         sock.precision(8);
         newly_opened = true;
      }
      sock << "mesh\n";

      mesh.Print(sock);

      if (newly_opened)
      {
         sock << "window_title '" << title << "'\n"
              << "window_geometry "
              << x << " " << y << " " << w << " " << h << "\n";
         if ( keys ) { sock << "keys " << keys << "\n"; }
         sock << endl;
      }

      connection_failed = !sock && !newly_opened;
   }
   while (connection_failed);
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

void VisualizeMesh(socketstream &sock, const char *vishost, int visport,
   ParMesh &pmesh, const char *title,
   int x, int y, int w, int h, const char *keys)
{
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
         sock << "mesh\n";
      }

      pmesh.PrintAsOne(sock);

      if (myid == 0 && newly_opened)
      {
         sock << "window_title '" << title << "'\n"
         << "window_geometry "
         << x << " " << y << " " << w << " " << h << "\n";
         if ( keys ) { sock << "keys " << keys << "\n"; }
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