#include "mfem.hpp"

using namespace mfem;

void VisualizeMesh(socketstream &sock, const char *vishost, int visport,
                   Mesh &mesh, const char *title,
                   int x = 0, int y = 0, int w = 400, int h = 400,
                   const char *keys = NULL);

int main(int argc, char *argv[])
{
   std::string mesh_name = "";
   bool viz = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_name, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&viz, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   args.ParseCheck();

   Mesh mesh(mesh_name.c_str());
   double hmin, hmax, kappamin, kappamax;
   Vector h_elem;
   mesh.GetCharacteristics(hmin, hmax, kappamin, kappamax, &h_elem);

   mfem::out << "Mesh: " << mesh_name << "\n"
             << "  Number of elements: " << mesh.GetNE() << "\n"
             << "  Number of vertices: " << mesh.GetNV() << "\n"
             << "  Min element size: " << hmin << "\n"
             << "  Max element size: " << hmax << "\n"
             << "  Avg element size: " << h_elem.Sum() / h_elem.Size() << "\n";
            //  << "  Min curvature radius: " << kappamin << "\n"
            //  << "  Max curvature radius: " << kappamax << "\n";

   if (viz)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sout(vishost, visport);
      VisualizeMesh(sout,vishost, visport, mesh, mesh_name.c_str(), 0, 0, 500, 500,"RjlmAb***********e");
   }

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
         sock << std::endl;
      }

      connection_failed = !sock && !newly_opened;
   }
   while (connection_failed);
}