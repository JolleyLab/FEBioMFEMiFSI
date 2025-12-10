#include "mfem.hpp"

using namespace mfem;

#ifdef MFEM_USE_MPI
void VisualizeMesh(socketstream &sock, const char *vishost, int visport,
                   ParMesh &pmesh, const char *title,
                   int x = 0, int y = 0, int w = 400, int h = 400,
                   const char *keys = NULL);
#else
void VisualizeMesh(socketstream &sock, const char *vishost, int visport,
                   Mesh &mesh, const char *title,
                   int x = 0, int y = 0, int w = 400, int h = 400,
                   const char *keys = NULL);
#endif

int main(int argc, char *argv[])
{
   #ifdef MFEM_USE_MPI
   Mpi::Init(argc, argv);
   Hypre::Init();
   #endif

   std::string mesh_name = "";
   bool vis = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_name, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&vis, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   args.ParseCheck();
   
   Mesh mesh(mesh_name.c_str());
   #ifdef MFEM_USE_MPI
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();
   #endif
   
   double hmin, hmax, havg, kappamin, kappamax;
   hmin = hmax = havg = kappamin = kappamax = 0.0;
   int num_elem, num_vert;
   Vector h_elem;
   #ifdef MFEM_USE_MPI   
   // num elem
   num_elem = pmesh.GetNE();
   MPI_Allreduce(MPI_IN_PLACE, &num_elem, 1, MPI_INT, MPI_SUM, pmesh.GetComm());

   // num nodes
   pmesh.EnsureNodes();
   auto nodes_space = pmesh.GetNodes()->FESpace();   
   auto par_nodes_space = dynamic_cast<ParFiniteElementSpace*>(nodes_space);
   int vdim = nodes_space->GetVDim();
   num_vert = par_nodes_space->GlobalTrueVSize() / vdim;

   // mesh sizes
   pmesh.GetCharacteristics(hmin, hmax, kappamin, kappamax);
   // compute havg
   double h_sum = 0.0;
   int dim = pmesh.Dimension();
   int sdim = pmesh.SpaceDimension();
   DenseMatrix J(sdim,dim);
   for (int i = 0; i < pmesh.GetNE(); i++)
   {
      pmesh.GetElementJacobian(i,J);
      h_sum += pow(fabs(J.Weight()), 1.0/real_t(dim));
   }
   MPI_Allreduce(MPI_IN_PLACE, &h_sum, 1, MPI_DOUBLE, MPI_SUM, pmesh.GetComm());
   havg = h_sum / num_elem;   
   #else
   mesh.GetCharacteristics(hmin, hmax, kappamin, kappamax, &h_elem);
   num_elem = mesh.GetNE();
   num_vert = mesh.GetNV();
   havg = h_elem.Sum() / h_elem.Size();
   #endif
   
   #ifdef MFEM_USE_MPI
   if (Mpi::Root())
   #endif
   {
      mfem::out << "Mesh: " << mesh_name << "\n"
               << "  Number of elements: " << num_elem << "\n"
               << "  Number of vertices: " << num_vert << "\n"
               << "  Min element size: " << hmin << "\n"
               << "  Max element size: " << hmax << "\n"
               << "  Avg element size: " << havg << "\n";
               //  << "  Min curvature radius: " << kappamin << "\n"
               //  << "  Max curvature radius: " << kappamax << "\n";
   }

   if (vis)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sout;
      #ifdef MFEM_USE_MPI
      VisualizeMesh(sout, vishost, visport, pmesh, mesh_name.c_str(), 0, 0, 500, 500);
      #else
      VisualizeMesh(sout,vishost, visport, mesh, mesh_name.c_str(), 0, 0, 500, 500);
      #endif
   }

   return 0;
}

#ifdef MFEM_USE_MPI
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
#else
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
         sock << std::endl;
      }

      connection_failed = !sock && !newly_opened;
   }
   while (connection_failed);
}
#endif