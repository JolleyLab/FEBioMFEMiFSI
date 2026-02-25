// Immersed annular disk

#include "MFEMFSI/IncompressibleFluid/IncompressibleFluidSolver.hpp"
#include "MFEMFSI/ImmersedFSI/ImmersedSolid.hpp"
#include "MFEMFSI/ImmersedFSI/Materials/materials.hpp"

using namespace mfem;
using namespace ifs_solver;
using namespace std;

struct IFS_Settings
{
   double fluid_density = 1.0; // kg/m^3
   double dyn_viscosity = 1.0; // Pa*s
   double kinvis = dyn_viscosity / fluid_density;

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
   double atol = 1e-10;
   double rtol = 1e-6;
   int minIt = 2;
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
   
   // fluid mesh
   int dim = 2;
   int n = 64; // fixed for given solid mesh
} ctx;

struct ImmersedSolidSettings
{
   int order = ctx.order;

   // solid material properties
   double mu = 1.0; // Pa
   double solid_density = 1.0; // kg/m^3

   // penalty parameter
   // double kappa = mu/10.0; // penalty value (default value)
   double kappa = 0.0; // penalty value (no effect for this problem!)

   // int form_type = ImmersedSolid::form_type::TotalLagrangian;
   int form_type = ImmersedSolid::form_type::UpdatedLagrangian;
   int vol_type = NearlyIncompressibleMaterial::volumetric_energy::ST91;
} stx;

struct iFSI_Settings
{
   double dt = 1e-3;
   int step_save = 1;

   int min_it = 2;
   int max_it = 10;

   double t_final = 10*dt; // 10 time steps

   bool visualization = false;
   int orderingtype = Ordering::byVDIM;
   // int orderingtype = Ordering::byNODES;

   // parameters
   double R = 0.25;
   double w = 0.06250;
   double len = 1.0;

   bool subtract_visc_fluid = false; // needed for this problem to work!
} itx;

void vel_exact(const Vector &x, double t, Vector &u)
{
   u.SetSize(x.Size());
   u = 0.0;
}

double vel_zero(const Vector &x, double t)
{
   return 0.0;
}

void grad_vel_exact(const Vector &x, double t, DenseMatrix &gradu)
{
   gradu.SetSize(x.Size());
   gradu = 0.0;
}

double pres_exact(const Vector &x, double t)
{
   double p_return = 0.0;

   // compute r
   Vector x_shifted(x.Size());
   x_shifted = x;
   x_shifted(0) -= 0.5;
   x_shifted(1) -= 0.5;

   // compute magnitude of vector
   double r = x_shifted.Norml2();

   // compute pres
   if (r >= itx.R + itx.w)
   {
      p_return = - M_PI * stx.mu * ( (itx.R + itx.w) * (itx.R + itx.w) - itx.R * itx.R)/(2.0 * itx.len*itx.len);
   }
   else if (r > itx.R && r < (itx.R + itx.w))
   {
      p_return = stx.mu * log((itx.R + itx.w)/r) - M_PI * stx.mu * ( (itx.R + itx.w) * (itx.R + itx.w) - itx.R * itx.R)/(2.0 * itx.len*itx.len);
   }
   else if (r <= itx.R)
   {
      p_return = stx.mu * log(1.0 + itx.w/itx.R) - M_PI * stx.mu * ( (itx.R + itx.w) * (itx.R + itx.w) - itx.R * itx.R)/(2.0 * itx.len*itx.len);
   }
   else
   {
      MFEM_ABORT("Error: r is out of range");
   }

   return p_return;
}

double ComputeVectorH1Error(GridFunction *u_gf, VectorCoefficient *ex_u, MatrixCoefficient *ex_gradu, const IntegrationRule *irs[]);         

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

   // Device configuration
   Device device(ctx.device_config.c_str());
   if (Mpi::Root()) 
   {
      mfem::out << "Using device: " << ctx.device_config << std::endl;
      device.Print();
      mfem::out << std::endl;
   }

   //////////// Initialize and parameters ////////////////

   // Read the mesh from the given mesh file.
   Mesh *fluid_mesh = nullptr;
   Mesh *solid_mesh = nullptr;
   if (ctx.dim == 2)
   {
      fluid_mesh = new Mesh(Mesh::MakeCartesian2D(ctx.n, ctx.n, Element::QUADRILATERAL, false, itx.len, itx.len));
      std::string solid_mesh_name = "annular_solid_mesh_" + std::to_string(ctx.n) + ".msh";
      solid_mesh = new Mesh(solid_mesh_name, 1, 1);
   }
   else
   {
      double dz;
      if (ctx.n == 64) 
      { 
         dz = 0.015625;
      }
      else if (ctx.n == 16) 
      { 
         dz = 0.0625;
      } 
      else MFEM_ABORT("n must be 16 or 64 for 3D");

      fluid_mesh = new Mesh(Mesh::MakeCartesian3D(ctx.n, ctx.n, 1.0, Element::HEXAHEDRON, itx.len, itx.len, dz));
      std::string solid_mesh_name = "annular_solid_mesh_" + std::to_string(ctx.n) + "_3D.msh"; 
      solid_mesh = new Mesh(solid_mesh_name, 1, 1);
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
   IncompressibleFluidSolver *flowsolver = new IncompressibleFluidSolver(fluid_pmesh, ctx.order, ctx.equal_order ? ctx.order : ctx.order - 1, ctx.kinvis, ctx.fluid_density, device, ctx.orderingtype);

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

   // all homogenous dirichlet conditions for fluid velocity
   Array<int> bdr_attr(flowsolver->GetMesh()->bdr_attributes.Max());
   bdr_attr = 1;
   if (ctx.dim == 3)
   {
      bdr_attr[0] = 0;
      bdr_attr[5] = 0;
   }
   flowsolver->AddNoSlipBC(bdr_attr);

   // since all dirichlet, \int p d\Omega = 0 is enforced
   flowsolver->AddZeroMeanPressureCondition();

   if (ctx.dim == 3)
   {
      // apply no z velocity
      Array<int> z_attr(flowsolver->GetMesh()->bdr_attributes.Max());
      z_attr = 0;
      z_attr[0] = 1;
      z_attr[5] = 1;
      flowsolver->AddVelDirichletBC(vel_zero,z_attr,3); // vz = 0, vx,vy free
   }

   // solid material
   ConstantCoefficient init_rhos(stx.solid_density);
   // material coeffs
   ConstantCoefficient mu_coeff(stx.mu);
   ConstantCoefficient kappa_coeff(stx.kappa);
   
   // add immersed solid
   MaterialType *mat_type = new NearlyIncompressibleAnnularDisk(stx.form_type,solid_pmesh->Dimension(),stx.vol_type,kappa_coeff,mu_coeff,itx.w);
   flowsolver->AddImmersedSolid(solid_pmesh,stx.order,&init_rhos,mat_type,stx.form_type,itx.subtract_visc_fluid);

   ///////////////////////////////////////////////////////////////////////////////////
   // setup solvers
   flowsolver->Setup(itx.dt);

   ///////////////////////////////////////////////////////////////////////////////////
   // Setup visualization
   std::string folder_name;
   folder_name = "results_annular_disk_" + std::to_string(ctx.n);
   if (ctx.dim == 3) folder_name += "_3D";
   folder_name += "_np" + std::to_string(Mpi::WorldSize());

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
   fdc->SetCycle(0);
   fdc->SetTime(0.0);
   fdc->Save();

   // visualization
   ParGridFunction *d_gf = flowsolver->GetImmersedSolidDisplacement();
   ParGridFunction *vib_gf = flowsolver->GetImmersedSolidVelocity();
   // copy solid mesh for outputting on deformed state
   ParMesh current_solid_pmesh(*solid_pmesh);
   GridFunction temp_gf(d_gf->FESpace());
   temp_gf = 0.0;

   std::string solid_filename;
   solid_filename = folder_name + "/solid";
   ParaViewDataCollection *sdc = new ParaViewDataCollection(solid_filename,&current_solid_pmesh);  
   sdc->SetDataFormat(VTKFormat::BINARY);
   sdc->SetHighOrderOutput(true);
   sdc->SetLevelsOfDetail(stx.order);
   sdc->SetPrecision(8);
   sdc->RegisterField("Displacement", d_gf);
   sdc->RegisterField("Velocity", vib_gf);
   sdc->SetCycle(0);
   sdc->SetTime(0.0);
   sdc->Save();

   ///////////////////////////////////////////////////////////////////////////////////
   socketstream sout, sout2, sout3;
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
         VisualizeField(sout3, vishost, visport, *d_gf, "Displacement",
            Wx+1000, Wy, Ww, Wh, "Rcj");
         // VisualizeMesh(sout3, vishost, visport, current_solid_pmesh, "Solid Mesh",
         //    Wx+1000, Wy, Ww, Wh, "pppppppma");
      }
      else
      {
         VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                  Wx, Wy, Ww, Wh, "Rcppppppppppppppppppppppppppppppppp");           
         VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                     Wx+500, Wy, Ww, Wh, "Rcppppppppppppppppppppppppppppppppp");
         VisualizeField(sout3, vishost, visport, *d_gf, "Displacement",
                        Wx+1000, Wy, Ww, Wh, "Rcpppppppppppppppppppppppppppppppppm");
         // VisualizeMesh(sout3, vishost, visport, current_solid_pmesh, "Solid Mesh",
         // Wx+1000, Wy, Ww, Wh, "Rmmma");        
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
         fdc->SetCycle(step);
         fdc->SetTime(t);
         fdc->Save();

         if (sdc)
         {
            // deform current solid mesh for visualization
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
            VisualizeField(sout3, vishost, visport, *d_gf, "Displacement",
                        Wx+1000, Wy, Ww, Wh, "Rcj");
            // VisualizeMesh(sout3, vishost, visport, current_solid_pmesh, "Solid Mesh",
            //    Wx+1000, Wy, Ww, Wh, "pppppppma");
         }
         else
         {
            VisualizeField(sout, vishost, visport, *u_gf, "Velocity",
                     Wx, Wy, Ww, Wh, "Rcppppppppppppppppppppppppppppppppp");           
            VisualizeField(sout2, vishost, visport, *p_gf, "Pressure",
                        Wx+500, Wy, Ww, Wh, "Rcppppppppppppppppppppppppppppppppp");
            VisualizeField(sout3, vishost, visport, *d_gf, "Displacement",
               Wx+1000, Wy, Ww, Wh, "Rcpppppppppppppppppppppppppppppppppm");
            // VisualizeMesh(sout3, vishost, visport, current_solid_pmesh, "Solid Mesh",
            // Wx+1000, Wy, Ww, Wh, "Rmmma");  
         } 
         MPI_Barrier(fluid_pmesh->GetComm()); 
      }
   }

   // print out relative L2 error at final time step
   double err_u = 0.0;
   double err_p = 0.0;
   VectorFunctionCoefficient vel_excoeff(fluid_pmesh->Dimension(), vel_exact);
   FunctionCoefficient p_excoeff(pres_exact);
   vel_excoeff.SetTime(t);
   p_excoeff.SetTime(t);
   
   // compute relative L2 error on finer mesh
   int order_err = ctx.order + 3;
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i=0; i < Geometry::NumGeom; ++i)
   {
      irs[i] = &(IntRules.Get(i, order_err));
   }
   err_u  = u_gf->ComputeL2Error(vel_excoeff, irs);
   err_p  = p_gf->ComputeL2Error(p_excoeff, irs);

   // compute velocity h1 norm
   MatrixFunctionCoefficient gradu_excoeff(fluid_pmesh->Dimension(),grad_vel_exact);

   // set time on coefficient
   gradu_excoeff.SetTime(t);
   double h1_err_u = ComputeVectorH1Error(u_gf,&vel_excoeff,&gradu_excoeff,irs);

   if (Mpi::Root())
   {
      std::cout << setprecision(16);
      // absolute h1 error
      // h1 norm for velocity
      std::cout << "|| u_h - u_ex ||_h1 = " << h1_err_u << "\n";

      // absolute L2 error
      std::cout << "|| u_h - u_ex || = " << err_u << "\n";
      std::cout << "|| p_h - p_ex || = " << err_p << "\n";  
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