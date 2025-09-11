#include "MFEMParExport.h"

// FEBio Headers
#include <XML/XMLWriter.h>
#include <FEBioLib/FEBioModel.h>
#include <FECore/FEPlotDataStore.h>
#include <FECore/fecore_enum.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEInitialCondition.h>
#include <FECore/log.h>
#include <XML/XMLReader.h>

// MFEM Headers
#include <mfem.hpp>
#include "MFEMImmersedElasticDomain.h"
#include "MFEMBoundaryCondition.h"
#include "MFEMLoad.h"

#include <fstream>
#include <sstream>

MFEMParExport::MFEMParExport()
{
	// Constructor implementation
}

MFEMParExport::~MFEMParExport()
{
   // Destructor implementation
}

void MFEMParExport::WriteFiles(std::string &feb_filename_, std::string &part_dir_, int np, bool immersed_febio_solid_)
{
   myid = 0;
   feb_filename = feb_filename_;
   part_dir = part_dir_;
   immersed_febio_solid = immersed_febio_solid_;

   // initialize reader to open original feb file
   std::string file_name = feb_filename + ".feb";
   bool file_opened = xml_reader.Open(file_name.c_str());
   MFEM_ASSERT(file_opened, "Error: could not open " + file_name);

   WriteMainFebFile();

   FEModel* fem = GetFEModel();
   FEMesh& femesh = fem->GetMesh();

   for (int i=0; i<np; ++i)
   {
      myid = i;
      WriteLocalFebFiles();
      
      // Reactivate all domains for next partition
      for (int d=0; d < femesh.Domains(); ++d)
      {
         FEDomain& dom = femesh.Domain(d);
         MFEMImmersedElasticDomain* im_dom = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
         if (im_dom)
         {
            im_dom->SetActive(true);
         }
      }
   }

   // close reader
   xml_reader.Close();
}

void MFEMParExport::WriteMainFebFile()
{
   std::string main_feb_file = part_dir + feb_filename + "_Main.feb";
   xml_writer.open(main_feb_file.c_str());
   
   // set modified formats
   XMLElement::intFormat = "%d";

   XMLElement el;

   // output root element
   el.name("febio_spec");
   el.add_attribute("version", "4.0");

   xml_writer.add_branch(el);
   {
      // module
      WriteModuleSection();

      // control
      WriteControlSection();

      // globals
      WriteGlobalsSection();

      // materials
      WriteMaterialsSection();

      // output
      WriteOutputSection();
   }
   xml_writer.close_branch(); // febio_spec

   // close the file
   xml_writer.close();
}

void MFEMParExport::WriteModuleSection()
{
   FEModel *fem = GetFEModel();
   XMLElement el;
   el.name("Module");
   el.add_attribute("type", fem->GetModuleName().c_str());
   xml_writer.add_branch(el);
   {
      // add module-specific elements here
      std::string units_str = fem->GetUnits() ? fem->GetUnits() : "";
      if (!units_str.empty())
      {
         auto write_units = [&](const std::string& u){ xml_writer.add_leaf("units", u.c_str()); };

         if (units_str == "2" || units_str == "SI")                  write_units("SI");
         else if (units_str == "3" || units_str == "mm-N-s")         write_units("mm-N-s");
         else if (units_str == "4" || units_str == "mm-kg-s")        write_units("mm-kg-s");
         else if (units_str == "5" || units_str == "um-nN-s")        write_units("um-nN-s");
         else if (units_str == "6" || units_str == "CGS")            write_units("CGS");
         else if (units_str == "7" || units_str == "mm-g-s")         write_units("mm-g-s");
         else if (units_str == "8" || units_str == "mm-mg-s")        write_units("mm-mg-s");
         else
         {
            // Fallback: write whatever string was provided
            xml_writer.add_leaf("units", units_str.c_str());
         }
      }
   }
   xml_writer.close_branch(); // module
}

void MFEMParExport::WriteControlSection()
{
   FEModel *fem = GetFEModel();
   int num_steps = fem->Steps();
   MFEM_ASSERT(num_steps == 1, "Error: TODO for more than one step");
   FEAnalysis *analysis = fem->GetCurrentStep();
   FESolver *solver = analysis->GetFESolver();
   xml_writer.add_branch("Control");
   {
      // Control parameters
      WriteParameterList(*analysis);

      // Solver
      XMLElement eli;
      eli.name("solver");
      solver->SetName("solver");
      WriteSectionComponent("Control", eli, *solver);
   }
   xml_writer.close_branch(); // Control
}

void MFEMParExport::WriteGlobalsSection()
{
   FEModel *fem = GetFEModel();
   XMLElement el;
   el.name("Globals");
   xml_writer.add_branch(el);
   {
      // add global variables here

      // Constants
      // create a list of strings for constant names
      // TODO: add more if needed, none are used by this plugin
      std::vector<std::string> constant_names = {"T","P","R","Fc"};
      xml_writer.add_branch("Constants");
      {
         for (const auto& name : constant_names)
         {
            double value = fem->GetGlobalConstant(name);
            xml_writer.add_leaf(name.c_str(), value);
         }
      }
      xml_writer.close_branch(); // Constants
   }
   xml_writer.close_branch(); // globals
}

void MFEMParExport::WriteMaterialsSection()
{
   FEModel *fem = GetFEModel();
   XMLElement el;
   el.name("Material");
   xml_writer.add_branch(el);
   {
      // add materials
      int num_mats = fem->Materials();
      for (int i = 0; i < num_mats; ++i)
      {
         FEMaterial* mat = fem->GetMaterial(i);
         {
            XMLElement eli;
            eli.name("material");
            WriteSectionComponent("Material", eli, *mat);
         }
      }
   }
   xml_writer.close_branch(); // Material

}

void MFEMParExport::WriteOutputSection()
{
   FEModel *fem = GetFEModel();
   XMLElement el;
   el.name("Output");
   xml_writer.add_branch(el);
   {
      FEPlotDataStore &plt_data_store = fem->GetPlotDataStore();
      XMLElement plt_el;
      plt_el.name("plotfile");
      plt_el.add_attribute("type", plt_data_store.GetPlotFileType().c_str());
      xml_writer.add_branch(plt_el);
      {
         int num_vars = plt_data_store.PlotVariables();
         for (int i=0; i < num_vars; i++)
         {
            const auto& var = plt_data_store.GetPlotVariable(i);
            XMLElement var_el;
            var_el.name("var");
            var_el.add_attribute("type", var.m_svar.c_str());
            xml_writer.add_empty(var_el);
         }
      }
      xml_writer.close_branch(); // plotfile
   }
   xml_writer.close_branch(); // Output
}

void MFEMParExport::WriteLocalFebFiles()
{
   std::ostringstream oss;
   oss << std::setw(10) << std::setfill('0') << myid;

   std::string local_feb_file = part_dir + feb_filename + "_Local_" + oss.str() + ".feb";
   xml_writer.open(local_feb_file.c_str());

   // set modified formats
   XMLElement::intFormat = "%d";

   XMLElement el;

   // output root element
   el.name("febio_spec");
   el.add_attribute("version", "4.0");

   xml_writer.add_branch(el);
   {
      // initial conditions
      WriteInitialConditionsSection();

      // boundary conditions
      WriteBoundaryConditionsSection();

      // loads
      WriteLoadsSection();

      if (immersed_febio_solid)
      {
         // mesh
         WriteMeshSection();

         // mesh domain
         WriteMeshDomainsSection();
      }
   }
   xml_writer.close_branch(); // febio_spec

   // close the file
   xml_writer.close();
}

void MFEMParExport::WriteInitialConditionsSection()
{
   FEModel *fem = GetFEModel();
   XMLElement el;
   el.name("Initial");
   xml_writer.add_branch(el);
   {
      // add initial conditions
      int num_ICs = fem->InitialConditions();
      for (int i = 0; i < num_ICs; ++i)
      {
         FEInitialCondition* ic = fem->InitialCondition(i);
         XMLElement ic_el("ic");
         WriteSectionComponent("Initial", ic_el, *ic);
      }
   }
   xml_writer.close_branch(); // Initial
}

void MFEMParExport::WriteBoundaryConditionsSection()
{
   FEModel *fem = GetFEModel();
   XMLElement el;
   el.name("Boundary");
   xml_writer.add_branch(el);
   {
      // add boundary conditions
      int num_BCs = fem->BoundaryConditions();
      for (int i = 0; i < num_BCs; ++i)
      {
         MFEMBoundaryCondition* bc = dynamic_cast<MFEMBoundaryCondition*>(fem->BoundaryCondition(i));
         MFEM_ASSERT(bc, "Error: MFEMBoundaryCondition not found");
         XMLElement bc_el("bc");
         bc_el.add_attribute("mfem_bdr_attr", bc->GetBoundaryID());
         WriteSectionComponent("Boundary", bc_el, *bc);
      }
   }
   xml_writer.close_branch(); // Boundary
}

void MFEMParExport::WriteLoadsSection()
{
   FEModel *fem = GetFEModel();
   XMLElement el;
   el.name("Loads");
   xml_writer.add_branch(el);
   {
      int num_loads = fem->ModelLoads();
      for (int i = 0; i < num_loads; ++i)
      {
         XMLElement load_el;

         // nodal loads
         MFEMNodalLoad* nl = dynamic_cast<MFEMNodalLoad*>(fem->ModelLoad(i));
         if (nl)
         {
            load_el.name("nodal_load");
            WriteSectionComponent("Loads", load_el, *nl);
         }

         // surface loads
         MFEMSurfaceLoad* sl = dynamic_cast<MFEMSurfaceLoad*>(fem->ModelLoad(i));
         if (sl)
         {
            load_el.name("surface_load");
            load_el.add_attribute("mfem_bdr_attr", sl->GetBoundaryID());
            WriteSectionComponent("Loads", load_el, *sl);
         }

         // body loads
         MFEMBodyLoad* bl = dynamic_cast<MFEMBodyLoad*>(fem->ModelLoad(i));
         if (bl)
         {
            load_el.name("body_load");
            load_el.add_attribute("mfem_dom_attr", bl->GetDomainAttr());
            WriteSectionComponent("Loads", load_el, *bl);
         }      
         
         // immersed surface loads
         FEImmersedSurfaceLoad* im_sl = dynamic_cast<FEImmersedSurfaceLoad*>(fem->ModelLoad(i));
         if (im_sl)
         {
            load_el.name("surface_load");
            load_el.add_attribute("mfem_bdr_attr", im_sl->GetBoundaryID());
            WriteSectionComponent("Loads", load_el, *im_sl);
         }
      }
   }
   xml_writer.close_branch(); // Loads
}

void MFEMParExport::WriteMeshSection()
{
   xml_writer.add_branch("Mesh");
   {
      WriteNodes();

      WriteElements();
   }
   xml_writer.close_branch(); // Mesh
}

// void MFEMParExport::WriteNodes()
// {
//    FEModel* fem = GetFEModel();
//    FEMesh& mesh = fem->GetMesh();

//    XMLElement tagNodes("Nodes");
//    xml_writer.add_branch(tagNodes);
//    {
//       XMLElement el("node");
//       int nid = el.add_attribute("id", 0);
//       for (int j = 0; j < mesh.Nodes(); ++j)
//       {
//          FENode& node = mesh.Node(j);
//          int m_nid = node.GetID();
//          el.set_attribute(nid, m_nid);
//          const vec3d r = node.m_r0;
//          std::ostringstream ss;
//          ss.setf(std::ios::scientific);
//          ss << std::setprecision(16)
//             << r.x << ", " << r.y << ", " << r.z;
//          el.value(ss.str().c_str());
//          xml_writer.add_leaf(el, false);
//       }
//    }
//    xml_writer.close_branch();
// }

// void MFEMParExport::WriteElements()
// {
//    FEModel* fem = GetFEModel();
//    FEMesh& mesh = fem->GetMesh();

//    int nset = 0;
//    int ncount = 0;
//    int nn[FEElement::MAX_NODES];
//    char szname[128] = { 0 };

//    // Assumes one combined immersed solid domain
//    bool immersed_solid_domain = false;
//    for (int d=0; d < mesh.Domains(); ++d)
//    {
//       FEDomain& dom = mesh.Domain(d);
//       MFEMImmersedElasticDomain* im_dom = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
//       if (im_dom && im_dom->IsActive())
//       {
//          immersed_solid_domain = true;

//          XMLElement xe("Elements");
//          xe.add_attribute("name", im_dom->GetName().c_str());
//          xml_writer.add_branch(xe);
//          {
//             for (int i=0; i<im_dom->Elements(); ++i)
//             {
//                FEElement& el = im_dom->ElementRef(i);            
//                XMLElement xej("elem");
//                xej.add_attribute("id", el.GetID());
//                int ne = el.Nodes();
//                for (int k = 0; k < ne; ++k) nn[k] = mesh.Node(el.m_node[k]).GetID();
//                xej.value(nn, ne);
//                xml_writer.add_leaf(xej, false);
//             }
//             nset++;
//          }
//          xml_writer.close_branch(); // Elements
//          break;
//       }
//    }
//    MFEM_ASSERT(immersed_solid_domain,"No immersed solid domain found!");
// }

mfem::Mesh* MFEMParExport::LoadMeshPart(const std::string &prefix)
{
   std::ifstream f(mfem::MakeParFilename(prefix + "_", myid, ".mesh", 10));
   return new mfem::Mesh(f, 1, 0);
}

void MFEMParExport::WriteNodes()
{
   // read in parallel solid mesh for current rank
   mfem::Mesh *mesh = LoadMeshPart(part_dir + "solid");
   MFEM_ASSERT(mesh, "Error: failed to load parallel solid mesh for rank " << myid);
   // here we assume a linear mesh

   XMLElement tagNodes("Nodes");
   xml_writer.add_branch(tagNodes);
   {
      XMLElement el("node");
      int nid = el.add_attribute("id", 0);
      for (int j = 0; j < mesh->GetNV(); ++j)
      {
         el.set_attribute(nid, j);
         const double* vertices = mesh->GetVertex(j);
         const vec3d r = vec3d(vertices[0], vertices[1], vertices[2]);
         std::ostringstream ss;
         ss.setf(std::ios::scientific);
         ss << std::setprecision(16)
            << r.x << ", " << r.y << ", " << r.z;
         el.value(ss.str().c_str());
         xml_writer.add_leaf(el, false);
      }
   }
   xml_writer.close_branch();

   delete mesh;
}

void MFEMParExport::WriteElements()
{
   FEModel* fem = GetFEModel();
   FEMesh& femesh = fem->GetMesh();

   // read in parallel solid mesh for current rank
   mfem::Mesh *mesh = LoadMeshPart(part_dir + "solid");
   MFEM_ASSERT(mesh, "Error: failed to load parallel solid mesh for rank " << myid);

   int nset = 0;
   int ncount = 0;
   int nn[FEElement::MAX_NODES];
   char szname[128] = { 0 };

   // Assumes one combined immersed solid domain
   bool immersed_solid_domain = false;
   for (int d=0; d < femesh.Domains(); ++d)
   {
      FEDomain& dom = femesh.Domain(d);
      MFEMImmersedElasticDomain* im_dom = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (im_dom && im_dom->IsActive())
      {
         immersed_solid_domain = true;
         int dom_attr = dom_attr_map->operator[](im_dom);

         mfem::Array<int> mesh_attrs = mesh->attributes;
         if (mesh_attrs.Find(dom_attr) == -1) 
         {
            // this mesh partition does not contain this domain
            im_dom->SetActive(false);
            continue; 
         }

         XMLElement xe("Elements");
         xe.add_attribute("name", im_dom->GetName().c_str());
         xml_writer.add_branch(xe);
         {
            for (int i=0; i<mesh->GetNE(); ++i)
            {
               auto element = mesh->GetElement(i);
               if (element->GetAttribute() != dom_attr) continue;

               XMLElement xej("elem");
               xej.add_attribute("id", i+1); // FEBio ID = MFEM ID + 1

               const int *node_ids = element->GetVertices();
               int ne = element->GetNVertices();
               for (int k = 0; k < ne; ++k) nn[k] = node_ids[k];
               xej.value(nn, ne);
               xml_writer.add_leaf(xej, false);
            }
            nset++;
         }
         xml_writer.close_branch(); // Elements
      }
   }
   MFEM_ASSERT(immersed_solid_domain,"No immersed solid domain found!");

   delete mesh;
}

const char* MFEMParExport::ElementTypeString(int ntype)
{
	const char* sztype = 0;
	switch (ntype)
	{
	case FE_Element_Type::FE_TET4G1 :
   case FE_Element_Type::FE_TET4G4 :
      sztype = "tet4"  ; break;

	case FE_Element_Type::FE_TET5G4  :	
      sztype = "tet5"  ; break;

	case FE_Element_Type::FE_TET10G1 :	
   case FE_Element_Type::FE_TET10G4 :
   case FE_Element_Type::FE_TET10G8 :
      sztype = "tet10" ; break;

	case FE_Element_Type::FE_PENTA6G6 :	
      sztype = "penta6"; break;

	case FE_Element_Type::FE_HEX8G1 :
   case FE_Element_Type::FE_HEX8G8 :
      sztype = "hex8"  ; break;

	case FE_Element_Type::FE_HEX27G27 : 
      sztype = "hex27" ; break;

   case FE_Element_Type::FE_PYRA5G8 :	
      sztype = "pyra5"; break;


	case FE_Element_Type::FE_QUAD4G4 :
   case FE_Element_Type::FE_QUAD4G16 :
      sztype = "quad4" ; break;

	case FE_Element_Type::FE_TRI3G1  :	
   case FE_Element_Type::FE_TRI3G3  :
   case FE_Element_Type::FE_TRI3G7  :
      sztype = "tri3"  ; break;

   case FE_Element_Type::FE_TRI6G3  :	
   case FE_Element_Type::FE_TRI6G4  :
   case FE_Element_Type::FE_TRI6G7  :
      sztype = "tri6"  ; break;

   case FE_Element_Type::FE_QUAD9G9 :	
      sztype = "quad9" ; break;

	default:
		assert(false);
	}
	return sztype;
}

// void MFEMParExport::WriteMeshDomainsSection()
// {
//    FEModel* fem = GetFEModel();
//    FEMesh& mesh = fem->GetMesh();

//    xml_writer.add_branch("MeshDomains");
//    {
//       for (int i = 0; i < mesh.Domains(); ++i)
//       {
//          FEDomain& dom = mesh.Domain(i);
//          MFEMImmersedElasticDomain* im_dom = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
//          assert(im_dom);
//          XMLElement el;
//          if (dom.Class() == FE_DOMAIN_SOLID)
//          {
//             el.name("SolidDomain");
//             el.add_attribute("name", im_dom->GetName().c_str());
//             el.add_attribute("mat", im_dom->GetMaterial()->GetName().c_str());
//             xml_writer.add_empty(el);
//          }
//          else
//          {
//             MFEM_ABORT("Invalid domain type");
//          }
//       }
//    }
//    xml_writer.close_branch(); // MeshDomains
// }

void MFEMParExport::WriteMeshDomainsSection()
{
   FEModel* fem = GetFEModel();
   FEMesh& mesh = fem->GetMesh();

   xml_writer.add_branch("MeshDomains");
   {
      for (int i = 0; i < mesh.Domains(); ++i)
      {
         FEDomain& dom = mesh.Domain(i);
         MFEMImmersedElasticDomain* im_dom = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
         if (im_dom && im_dom->IsActive())
         {   
            XMLElement el;
            if (dom.Class() == FE_DOMAIN_SOLID)
            {
               el.name("SolidDomain");
               el.add_attribute("name", im_dom->GetName().c_str());
               el.add_attribute("mat", im_dom->GetMaterial()->GetName().c_str());
               xml_writer.add_empty(el);
            }
            else
            {
               MFEM_ABORT("Invalid domain type");
            }
         }
      }
   }
   xml_writer.close_branch(); // MeshDomains
}

void MFEMParExport::WriteComponent(FECoreBase &base, XMLElement &el)
{   
   // get the type string    
   const char* sztype = base.GetTypeStr(); assert(sztype);
   if (sztype && sztype[0])
   {
      // only add the type attribute if the tag name is different from the type string
      // there is an unfortunate coincidence that the nodal_load has the same type string
      // as tag name. We need to make an exception for this. 
      if ((strcmp(sztype, el.name()) != 0) || (strcmp(sztype, "nodal_load") == 0))
         el.add_attribute("type", sztype);
   }

   // see if there are any attribute parameters
   // Get the parameter list from the base
   FEParameterList &param_list = base.GetParameterList();
   list<FEParam>::iterator it = param_list.first();

   // Iterate through the parameters and write them
   int NP = base.Parameters();
   for (int i = 0; i < NP; i++)
   {
      FEParam &p = *it++;
      if (p.GetFlags() & FEParamFlag::FE_PARAM_ATTRIBUTE)
      {
         switch (p.type())
         {
         case FEParamType::FE_PARAM_INT:
         {
            if (p.enums())
            {
               const char* v = p.enumKey();
               if (v == nullptr) {
                  MFEM_ABORT("Invalid key for parameter %s" << p.name());
                  v = "(invalid)";
               }
               el.add_attribute(p.name(), v);
            }
            else el.add_attribute(p.name(), p.value<int>());
         }
         break;
         case FEParamType::FE_PARAM_STRING:
         {
            std::string s = p.value<std::string>();
            el.add_attribute(p.name(), s);
         }
         break;
         default:
            MFEM_ABORT("Unsupported attribute parameter type");
         }
      }
   }

   if ((base.Parameters() == 0) && (base.Properties() == 0))
   {
      xml_writer.add_empty(el);
      return;
   }
   
   xml_writer.add_branch(el);
   {
      // write parameters
      if (base.Parameters())
         WriteParameterList(base);

      // write properties
      int NP = base.Properties();
      for (int i = 0; i < NP; i++)
      {
         FECoreBase* prop = base.GetProperty(i);
         if (prop)
         {
            XMLElement eli(prop->GetName().c_str());
            std::string name = prop->GetName();
            if(!name.empty())
            {
               eli.add_attribute("name", name.c_str());
            }

            WriteComponent(*prop, eli);
         }         
      }
   }
   xml_writer.close_branch();
}

void MFEMParExport::WriteParameterList(FECoreBase &base)
{
   // Get the parameter list from the base
   FEParameterList &param_list = base.GetParameterList();
   list<FEParam>::iterator it = param_list.first();

   // Iterate through the parameters and write them
   int NP = base.Parameters();
   for (int i = 0; i < NP; i++)
   {
      FEParam &p = *it++;
      WriteParam(p);
   }
}

void MFEMParExport::WriteParam(FEParam &p)
{
   XMLElement el;
   el.name(p.name());

   // write
   switch (p.type())
   {
      case FEParamType::FE_PARAM_INT:
      {
         auto* v = static_cast<int*>(p.data_ptr());

         // check for the pre_partition_opt and change form 0 (partitioning) to 1 (reading)
         if (std::strcmp(p.name(), "pre_partition_opt") == 0)
         {
            el.value(1);
         }
         else el.value(*v);
         break;
      }
      case FEParamType::FE_PARAM_DOUBLE:
      {
         auto* v = static_cast<double*>(p.data_ptr());
         el.value(*v);
         break;
      }
      case FEParamType::FE_PARAM_BOOL:
      {
         auto* v = static_cast<bool*>(p.data_ptr());
         el.value((*v ? 1 : 0));
         break;
      }
      case FEParamType::FE_PARAM_STRING:
      {
         auto* v = static_cast<std::string*>(p.data_ptr());
         el.value(v->c_str());
         break;
      }
      case FEParamType::FE_PARAM_VEC3D:
      {
         const auto* v = static_cast<const vec3d*>(p.data_ptr());
         std::ostringstream ss;
         ss.setf(std::ios::scientific);
         ss << std::setprecision(16)
            << v->x << ", " << v->y << ", " << v->z;
         el.value(ss.str().c_str());  // pass as C-string
         break;
      }
      case FEParamType::FE_PARAM_VEC2D:
      {
         const auto* v = static_cast<const vec2d*>(p.data_ptr());
         std::ostringstream ss;
         ss.setf(std::ios::scientific);
         ss << std::setprecision(16)
            << v->x() << ", " << v->y();
         el.value(ss.str().c_str());
         break;
      }
      case FEParamType::FE_PARAM_MAT3D:
      {
         const auto* v = static_cast<const mat3d*>(p.data_ptr());
         std::ostringstream ss;
         ss.setf(std::ios::scientific);
         ss << std::setprecision(16)
            << (*v)(0,0) << ", " << (*v)(0,1) << ", " << (*v)(0,2) << ", "
            << (*v)(1,0) << ", " << (*v)(1,1) << ", " << (*v)(1,2) << ", "
            << (*v)(2,0) << ", " << (*v)(2,1) << ", " << (*v)(2,2);
         el.value(ss.str().c_str());
         break;
      }
      case FEParamType::FE_PARAM_MAT3DS:
      {
         const auto* v = static_cast<const mat3ds*>(p.data_ptr());
         const double m00 = v->xx(), m11 = v->yy(), m22 = v->zz();
         const double m01 = v->xy(), m02 = v->xz(), m12 = v->yz();
         std::ostringstream ss;
         ss.setf(std::ios::scientific);
         ss << std::setprecision(16)
            << m00 << ", " << m11 << ", " << m22 << ", "
            << m01 << ", " << m12 << ", " << m02;
         el.value(ss.str().c_str());
         break;
      }
      case FEParamType::FE_PARAM_VEC3D_MAPPED:
      {
         feLogWarning("Parameter type vec3d mapped not fully supported");
         const auto* v = static_cast<const vec3d*>(p.data_ptr());
         
         // FEParamVec3& v = p.value<FEParamVec3>();
         // FEVec3dValuator* vec_valuator = v.valuator();
         
         std::ostringstream ss;
         ss.setf(std::ios::scientific);
         ss << std::setprecision(16)
            << v->x << ", " << v->y << ", " << v->z;
         el.value(ss.str().c_str());
         break;
      }
      case FEParamType::FE_PARAM_DOUBLE_MAPPED:
      {
         feLogWarning("Parameter type double mapped not fully supported");
         FEParamDouble& v = p.value<FEParamDouble>();
         FEScalarValuator* scalar_valuator = v.valuator();

         if (scalar_valuator->isConst())
         {
            double const_value = *scalar_valuator->constValue();
            el.value(const_value);
         }
         else
         {
            MFEM_ABORT("Error reading mapped double parameter");
         }
         break;
      }
      default:
      {
         feLogError("Unsupported parameter type in control section");
         break;
      }
   }
   xml_writer.add_leaf(el);
}

// void MFEMParExport::WriteSectionComponent(const char* szname, XMLElement &el, FECoreBase &base)
// {
//    // look for the tag
//    XMLTag tag;
//    std::string path_name = std::string("febio_spec/") + szname;
//    if (xml_reader.FindTag(path_name.c_str(), tag) == false)
//    {
//       MFEM_ABORT("Error: could not read " + path_name);
//    }

//    // make sure this tag has children
//    if (tag.isleaf())
//    {
//       // xml_writer.add_empty(el);
//       MFEM_ABORT("TODO: empty tag");
//       return;
//    } 

//    // search for the correct tag
//    ++tag;
//    bool found_tag = false;
//    std::ifstream *fs = xml_reader.GetFileStream();
//    do
//    {
//       const char* tag_name = tag.AttributeValue("name", true);
//       const char* sz_tag = tag.m_sztag.c_str();
//       const char* base_name = base.GetName().c_str();
//       if (tag_name && strcmp(tag_name, base_name) == 0)
//       {
//          found_tag = true;
//          break;
//       }
//       else if (sz_tag && strcmp(sz_tag, base_name) == 0)
//       {
//          found_tag = true;
//          break;
//       }
//       else 
//       {
//          xml_reader.SkipTag(tag);
//          ++tag;
//       }
//    }
//    while(found_tag == false && fs->eof() == false);
//    MFEM_ASSERT(found_tag, "Error: could not find tag " + std::string(el.name()));

//    ParseTag(tag, el, base);

//    xml_reader.Close();
// }

void MFEMParExport::WriteSectionComponent(const char* szname, XMLElement &el, FECoreBase &base)
{
   // Build cache once (skips large Mesh section)
   BuildTopLevelSectionCache();

   std::string section_xml;
   if (!GetCachedSection(szname, section_xml))
   {
      MFEM_ABORT("Section not cached (or Mesh) : " << szname);
      return;
   }

   // Wrap in artificial root so we can reuse current parser logic safely.
   std::string wrapped = "<root>\n" + section_xml + "\n</root>\n";

   // Write wrapped to a temporary in-memory file substitute:
   // Simplest approach: write to a std::stringstream then dump to a temp file.
   std::string tmpFile = part_dir + "_tmp_cached_section.xml";
   {
      std::ofstream ofs(tmpFile, std::ios::binary);
      ofs << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << wrapped;
   }

   XMLReader local;
   if (!local.Open(tmpFile.c_str()))
   {
      MFEM_ABORT("Failed to open temp cached section file");
      return;
   }

   XMLTag root;
   if (!local.FindTag("root", root))
   {
      local.Close();
      MFEM_ABORT("Temp root not found");
      return;
   }

   // Iterate children (top-level: only one target section tag)
   ++root;
   bool sectionFound = false;
   while (!root.isend())
   {
      if (strcmp(root.Name(), szname) == 0)
      {
         sectionFound = true;
         break;
      }
      local.SkipTag(root);
      ++root;
   }
   MFEM_ASSERT(sectionFound, "Cached section root mismatch");

   // Now search for component inside that section
   if (!root.isleaf())
   {
      ++root;
      bool found_tag = false;
      const char* base_name = base.GetName().c_str();
      // MFEM_ASSERT(base_name && base_name[0], "Error: component has no name");
      while (!root.isend())
      {
         const char* tag_attr_name = root.AttributeValue("name", true);
         const char* tag_label     = root.m_sztag.c_str();
         if ((tag_attr_name && strcmp(tag_attr_name, base_name) == 0))
         {
            found_tag = true;
            break;
         }
         else if (tag_label && strcmp(tag_label, base_name) == 0)
         {
            found_tag = true;
            break;
         }
         local.SkipTag(root);
         ++root;
      }
      MFEM_ASSERT(found_tag, std::string("Component not found in cached section: ") + base_name);

      // Parse into output XML
      ParseTag(root, el, base);
   }
   else
   {
      // Leaf section (unlikely for sections used here)
      xml_writer.add_empty(el);
      MFEM_ABORT("Error: empty tag");
   }

   local.Close();

   // remove temp file
   std::remove(tmpFile.c_str());
}

void MFEMParExport::BuildTopLevelSectionCache()
{
   if (m_sectionCacheBuilt) return;
   m_sectionCache.clear();

   const std::string file_name = feb_filename + ".feb";
   std::ifstream fs(file_name, std::ios::binary);
   if(!fs) return;

   // Read whole file once
   std::string file;
   fs.seekg(0, std::ios::end);
   file.resize(static_cast<size_t>(fs.tellg()));
   fs.seekg(0, std::ios::beg);
   fs.read(file.data(), static_cast<std::streamsize>(file.size()));

   // Find <febio_spec
   size_t root_pos = file.find("<febio_spec");
   if (root_pos == std::string::npos) { m_sectionCacheBuilt = true; return; }

   // Move to end of root start tag '>'
   size_t pos = file.find('>', root_pos);
   if (pos == std::string::npos) { m_sectionCacheBuilt = true; return; }
   ++pos;

   // Simple top-level scan: depth relative to inside <febio_spec> content
   while (pos < file.size())
   {
      // skip whitespace
      while (pos < file.size() && isspace(static_cast<unsigned char>(file[pos]))) ++pos;
      if (pos >= file.size()) break;

      if (file.compare(pos, 11, "</febio_spec") == 0) break;
      if (file[pos] != '<') { ++pos; continue; }
      if (file.compare(pos, 4, "<!--") == 0) {
         // skip comment
         size_t cend = file.find("-->", pos+4);
         if (cend == std::string::npos) break;
         pos = cend + 3;
         continue;
      }
      size_t tag_start = pos;
      ++pos;
      bool closing = false;
      if (file[pos] == '/') { closing = true; }

      if (closing) {
         // unexpected closing at top-level -> done
         size_t close_gt = file.find('>', pos);
         if (close_gt == std::string::npos) break;
         pos = close_gt + 1;
         continue;
      }

      // Parse tag name
      size_t name_begin = pos;
      while (pos < file.size() && !isspace(static_cast<unsigned char>(file[pos])) &&
             file[pos] != '>' && file[pos] != '/') ++pos;
      std::string section_name = file.substr(name_begin, pos - name_begin);

      // We are at start tag of a top-level child (<Section ...>)
      // Find end of start tag
      size_t start_tag_end = file.find('>', pos);
      if (start_tag_end == std::string::npos) break;
      bool self_closing = (file[start_tag_end-1] == '/');

      if (section_name == "Mesh")
      {
         // Skip entire Mesh subtree by depth count
         if (self_closing) { pos = start_tag_end + 1; continue; }
         size_t scan = start_tag_end + 1;
         int depth = 1;
         while (scan < file.size() && depth > 0)
         {
            size_t lt = file.find('<', scan);
            if (lt == std::string::npos) break;
            if (file.compare(lt, 4, "<!--") == 0) {
               size_t cend = file.find("-->", lt+4);
               if (cend == std::string::npos) break;
               scan = cend + 3;
               continue;
            }
            bool closing2 = (lt + 1 < file.size() && file[lt+1] == '/');
            size_t gt = file.find('>', lt+1);
            if (gt == std::string::npos) break;
            if (!closing2) {
               // check self closing
               if (file[gt-1] != '/') depth++;
            } else {
               depth--;
            }
            scan = gt + 1;
         }
         pos = scan;
         continue;
      }

      // Capture non-Mesh section
      size_t section_end_pos = start_tag_end + 1;
      if (!self_closing)
      {
         // depth count just for this section
         size_t scan = start_tag_end + 1;
         int depth = 1;
         while (scan < file.size() && depth > 0)
         {
            size_t lt = file.find('<', scan);
            if (lt == std::string::npos) break;
            if (file.compare(lt, 4, "<!--") == 0) {
               size_t cend = file.find("-->", lt+4);
               if (cend == std::string::npos) break;
               scan = cend + 3;
               continue;
            }
            bool closing2 = (lt + 1 < file.size() && file[lt+1] == '/');
            size_t gt = file.find('>', lt+1);
            if (gt == std::string::npos) break;
            if (!closing2) {
               if (file[gt-1] != '/') depth++;
            } else {
               depth--;
            }
            scan = gt + 1;
            if (depth == 0) section_end_pos = scan;
         }
         pos = section_end_pos;
      }
      else
      {
         pos = start_tag_end + 1;
      }

      // Store raw substring of this section
      m_sectionCache[section_name] = file.substr(tag_start, section_end_pos - tag_start);
   }

   m_sectionCacheBuilt = true;
}

bool MFEMParExport::GetCachedSection(const std::string& name, std::string& out) const
{
   auto it = m_sectionCache.find(name);
   if (it == m_sectionCache.end()) return false;
   out = it->second;
   return true;
}

void MFEMParExport::ParseTag(XMLTag &tag, XMLElement &el, FECoreBase &base)
{
   WriteAttributes(tag, el, base);

   // write the element
   xml_writer.add_branch(el,false);
   {
      XMLElement eli;

      // write parameters
      if (!tag.isleaf())
      {
         ++tag;
         do
         {
            if (tag.isleaf())
            {
               // write the leaf
               eli.name(tag.Name());
               WriteAttributes(tag,eli,base);
               WriteParameters(tag,eli,base);

               if (strcmp(eli.name(), "pre_partition_opt") == 0)
               {
                  // special case for pre_partition_opt
                  eli.value(1); // change from 0 (partitioning) to 1 (reading)
               }

               xml_writer.add_leaf(eli);
            }
            else if (tag.isempty())
            {
               // xml_writer.add_empty(eli);
               xml_writer.close_branch();
               return;
            }
            else
            {
               // recursion
               ParseTag(tag,eli,base);
            }
            ++tag;
         }
         while (!tag.isend());
      }
      else
      {
         // xml_writer.add_empty(eli);
         xml_writer.close_branch();
         return;
      }
   }
   xml_writer.close_branch();
}

void MFEMParExport::WriteAttributes(XMLTag &tag, XMLElement &el, FECoreBase &base)
{
   int Num_Attr = tag.m_att.size();
   for (int i = 0; i < Num_Attr; ++i)
   {
      XMLAtt& att = tag.m_att[i];
      const char* szatt = att.name();
      const char* szval = att.cvalue();
      if (szatt && szval)
      {
         el.add_attribute(szatt, szval);
      }
   }
}

void MFEMParExport::WriteParameters(XMLTag& tag, XMLElement &el, FECoreBase &base)
{
   el.value(tag.m_szval);
}