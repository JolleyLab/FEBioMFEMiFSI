#include "mfem.hpp"

using namespace mfem;

#ifdef MFEM_USE_MPI
void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    const ParGridFunction &gf, const char *title,
                    int x = 0, int y = 0, int w = 400, int h = 400,
                    const char *keys = NULL, bool vec = false);

void ComputeElementSizes(ParMesh &mesh, double &hmin, double &hmax, double &havg, int &num_elem, int &num_vert);
#else
void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    GridFunction &gf, const char *title,
                    int x = 0, int y = 0, int w = 400, int h = 400,
                    const char *keys = NULL, bool vec = false);

void ComputeElementSizes(Mesh &mesh, double &hmin, double &hmax, double &havg, int &num_elem, int &num_vert);
#endif

int main(int argc, char *argv[])
{
   #ifdef MFEM_USE_MPI
   Mpi::Init(argc, argv);
   Hypre::Init();
   #endif

   std::string fluid_mesh_name = "";
   std::string solid_mesh_name = "";
   bool vis = true;
   bool save_to_paraview = false;

   OptionsParser args(argc, argv);
   args.AddOption(&fluid_mesh_name, "-fm", "--fluid-mesh",
                  "Fluid mesh file to use.");
   args.AddOption(&solid_mesh_name, "-sm", "--solid-mesh",
                  "Solid mesh file to use.");
   args.AddOption(&save_to_paraview, "-sp", "--save-to-paraview", "-no-sp",
                  "--no-save-to-paraview",
                  "Save the mesh size ratio field to a Paraview file.");
   args.AddOption(&vis, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   args.ParseCheck();

   Mesh fluid_mesh(fluid_mesh_name.c_str());
   int dim = fluid_mesh.Dimension();
   int sdim = fluid_mesh.SpaceDimension();
   #ifdef MFEM_USE_MPI
   ParMesh fluid_pmesh(MPI_COMM_WORLD, fluid_mesh);
   fluid_mesh.Clear();
   #endif
   double fluid_hmin, fluid_hmax, fluid_havg;
   int fluid_num_elem, fluid_num_vert;
   #ifdef MFEM_USE_MPI
   ComputeElementSizes(fluid_pmesh, fluid_hmin, fluid_hmax, fluid_havg,
                       fluid_num_elem, fluid_num_vert);
   #else
   ComputeElementSizes(fluid_mesh, fluid_hmin, fluid_hmax, fluid_havg,
                       fluid_num_elem, fluid_num_vert);
   #endif
   
   #ifdef MFEM_USE_MPI
   if (Mpi::Root())
   #endif
   {
      mfem::out << "Mesh: " << fluid_mesh_name << "\n"
             << "  Number of elements: " << fluid_num_elem << "\n"
             << "  Number of vertices: " << fluid_num_vert << "\n"
             << "  Min element size: " << fluid_hmin << "\n"
             << "  Max element size: " << fluid_hmax << "\n"
             << "  Avg element size: " << fluid_havg << "\n";
   }   

   Mesh solid_mesh(solid_mesh_name.c_str());
   #ifdef MFEM_USE_MPI
   ParMesh solid_pmesh(MPI_COMM_WORLD, solid_mesh);
   solid_mesh.Clear();
   #endif
   double solid_hmin, solid_hmax, solid_havg;
   int solid_num_elem, solid_num_vert;
   #ifdef MFEM_USE_MPI
   ComputeElementSizes(solid_pmesh, solid_hmin, solid_hmax, solid_havg,
                       solid_num_elem, solid_num_vert);
   #else
   ComputeElementSizes(solid_mesh, solid_hmin, solid_hmax, solid_havg,
                       solid_num_elem, solid_num_vert);
   #endif
   #ifdef MFEM_USE_MPI
   if (Mpi::Root())
   #endif
   {
      mfem::out << "Mesh: " << solid_mesh_name << "\n"
               << "  Number of elements: " << solid_num_elem << "\n"
               << "  Number of vertices: " << solid_num_vert << "\n"
               << "  Min element size: " << solid_hmin << "\n"
               << "  Max element size: " << solid_hmax << "\n"
               << "  Avg element size: " << solid_havg << "\n";

      mfem::out << "Ratio of fluid element size to max solid element size:\n"
               << "  Min: " << fluid_hmin / solid_hmax << "\n"
               << "  Max: " << fluid_hmax / solid_hmax << "\n"
               << "  Avg: " << fluid_havg / solid_hmax << "\n";
   }
   
   // Visualize ratio of fluid element size to max solid element size
   L2_FECollection fec(0, dim);
   #ifdef MFEM_USE_MPI
   ParFiniteElementSpace fes(&fluid_pmesh, &fec, 1);
   ParGridFunction h_ratio(&fes);
   #else
   FiniteElementSpace fes(&fluid_mesh, &fec, 1);
   GridFunction h_ratio(&fes);
   #endif
   h_ratio = 0.0;
   DenseMatrix J(sdim, dim);
   double h_elem;
   Array<int> fluid_elem_dofs;
   #ifdef MFEM_USE_MPI
   for (int i = 0; i < fluid_pmesh.GetNE(); i++)
   {  
      fluid_pmesh.GetElementJacobian(i,J);
      h_elem = pow(fabs(J.Weight()), 1.0/real_t(dim));
      fes.GetElementDofs(i, fluid_elem_dofs);
      MFEM_ASSERT(fluid_elem_dofs.Size() == 1, "Unexpected dof size");
      h_ratio[fluid_elem_dofs[0]] = h_elem / solid_hmax;
   }
   #else
   for (int i = 0; i < fluid_mesh.GetNE(); i++)
   {  
      fluid_mesh.GetElementJacobian(i,J);
      h_elem = pow(fabs(J.Weight()), 1.0/real_t(dim));
      fes.GetElementDofs(i, fluid_elem_dofs);
      MFEM_ASSERT(fluid_elem_dofs.Size() == 1, "Unexpected dof size");
      h_ratio[fluid_elem_dofs[0]] = h_elem / solid_hmax;
   }
   #endif

   if (save_to_paraview)
   {
      ParaViewDataCollectionBase *pvdc = nullptr;
      #ifdef MFEM_USE_MPI
      #ifdef MFEM_USE_HDF5
      pvdc = new ParaViewHDFDataCollection("MeshSizeRatio", &fluid_pmesh);
      #else
      pvdc = new ParaViewDataCollection("MeshSizeRatio", &fluid_pmesh);
      #endif       
      #else
      pvdc = new ParaViewDataCollection("MeshSizeRatio", &fluid_mesh);
      #endif
      pvdc->SetDataFormat(mfem::VTKFormat::BINARY32);
      pvdc->SetHighOrderOutput(false);
      pvdc->SetLevelsOfDetail(1);
      pvdc->SetCompression(true);
      pvdc->SetCycle(0);
      pvdc->SetTime(0);
      pvdc->RegisterField("MeshSizeRatio", &h_ratio);

      pvdc->Save();
      delete pvdc;
   }

   if (vis)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sout;
      VisualizeField(sout, vishost, visport, h_ratio,
                     "Fluid to Solid Mesh Size Ratio",
                     0, 0, 600, 600,
                     "ci", false);
   }

   return 0;
}

#ifdef MFEM_USE_MPI
void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    const ParGridFunction &gf, const char *title,
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
         sock << std::endl;
      }

      if (myid == 0)
      {
         connection_failed = !sock && !newly_opened;
      }
      MPI_Bcast(&connection_failed, 1, MPI_INT, 0, comm);
   }
   while (connection_failed);
}

void ComputeElementSizes(ParMesh &mesh, double &hmin, double &hmax, double &havg, int &num_elem, int &num_vert)
{
   // num elem
   num_elem = mesh.GetNE();
   MPI_Allreduce(MPI_IN_PLACE, &num_elem, 1, MPI_INT, MPI_SUM, mesh.GetComm());

   // num nodes
   mesh.EnsureNodes();
   auto nodes_space = mesh.GetNodes()->FESpace();   
   auto par_nodes_space = dynamic_cast<ParFiniteElementSpace*>(nodes_space);
   int vdim = nodes_space->GetVDim();
   num_vert = par_nodes_space->GlobalTrueVSize() / vdim;

   // mesh sizes
   double h_sum = 0.0;
   double h = 0.0;
   hmin = infinity();
   hmax = -hmin;
   int dim = mesh.Dimension();
   int sdim = mesh.SpaceDimension();
   DenseMatrix J(sdim,dim);
   for (int i = 0; i < mesh.GetNE(); i++)
   {
      mesh.GetElementJacobian(i,J);
      h = pow(fabs(J.Weight()), 1.0/real_t(dim));
      h_sum += h;

      if (h < hmin) { hmin = h; }
      if (h > hmax) { hmax = h; }
   }
   MPI_Allreduce(MPI_IN_PLACE, &hmin, 1, MPI_DOUBLE, MPI_MIN, mesh.GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &hmax, 1, MPI_DOUBLE, MPI_MAX, mesh.GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &h_sum, 1, MPI_DOUBLE, MPI_SUM, mesh.GetComm());
   havg = h_sum / num_elem;
}

#else
void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    GridFunction &gf, const char *title,
                    int x, int y, int w, int h, const char * keys, bool vec)
{
   Mesh &mesh = *gf.FESpace()->GetMesh();

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
      sock << "solution\n";

      mesh.Print(sock);
      gf.Save(sock);

      if (newly_opened)
      {
         sock << "window_title '" << title << "'\n"
              << "window_geometry "
              << x << " " << y << " " << w << " " << h << "\n";
         if ( keys ) { sock << "keys " << keys << "\n"; }
         else { sock << "keys maaAc\n"; }
         if ( vec ) { sock << "vvv"; }
         sock << std::endl;
      }

      connection_failed = !sock && !newly_opened;
   }
   while (connection_failed);
}

void ComputeElementSizes(Mesh &mesh, double &hmin, double &hmax, double &havg, int &num_elem, int &num_vert)
{
   Vector h_elem;
   double kappamin, kappamax;
   mesh.GetCharacteristics(hmin, hmax, kappamin, kappamax, &h_elem);
   havg = h_elem.Sum() / h_elem.Size();
   num_elem = mesh.GetNE();
   num_vert = mesh.GetNV();
}
#endif