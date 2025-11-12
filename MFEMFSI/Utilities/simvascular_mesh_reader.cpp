#include "simvascular_mesh_reader.hpp"
#include "../../Software/mfem/mesh/mesh_readers.cpp"

using namespace mfem;
using namespace std;

#ifdef MFEM_USE_ZLIB
void SimvascularMeshReader::Load(std::istream &input, int generate_edges, int refine, bool fix_orientation)
{
   SVLoader(input, generate_edges);
   Finalize(refine, fix_orientation);
}

SimvascularMeshReader::SimvascularMeshReader(const std::string &filename, int generate_edges, int refine, bool fix_orientation)
   : mfem::Mesh()
{
   named_ifgzstream imesh(filename);
   if (!imesh)
   {
      // Abort with an error message.
      MFEM_ABORT("Mesh file not found: " << filename << '\n');
   }
   else
   {
      Load(imesh, generate_edges, refine, fix_orientation);
   }
}

void SimvascularMeshReader::SVLoader(std::istream &input, int generate_edges)
{
   int curved = 0, read_gf = 1;
   bool finalize_topo = true;

   if (!input)
   {
      MFEM_ABORT("Input stream is not open");
   }

   Clear();

   string mesh_type;
   input >> ws;
   getline(input, mesh_type);
   filter_dos(mesh_type);
   
   if (mesh_type.rfind("<VTKFile ") == 0 || mesh_type.rfind("<?xml") == 0)
   {
      SV_ReadXML_VTKMesh(input, curved, read_gf, finalize_topo, mesh_type);
   }
   else
   {
      MFEM_ABORT("Simvascular mesh is expected to be a vtk xml mesh!" << mesh_type);
   }   

   // at this point the following should be defined:
   //  1) Dim
   //  2) NumOfElements, elements
   //  3) NumOfBdrElements, boundary
   //  4) NumOfVertices, with allocated space in vertices
   //  5) curved
   //  5a) if curved == 0, vertices must be defined
   //  5b) if curved != 0 and read_gf != 0,
   //         'input' must point to a GridFunction
   //  5c) if curved != 0 and read_gf == 0,
   //         vertices and Nodes must be defined
   // optional:
   //  1) el_to_edge may be allocated (as in the case of P2 VTK meshes)
   //  2) ncmesh may be allocated

   // FinalizeTopology() will:
   // - assume that generate_edges is true
   // - assume that refine is false
   // - does not check the orientation of regular and boundary elements
   if (finalize_topo)
   {
      // don't generate any boundary elements, especially in parallel
      bool generate_bdr = false;

      FinalizeTopology(generate_bdr);
   }

   if (curved && read_gf)
   {
      MFEM_ABORT("Curved meshes with read_gf != 0 are not supported");
   }

   // Finalize(...) should be called after this, if needed.
}

void SimvascularMeshReader::SV_ReadXML_VTKMesh(std::istream &input, int &curved, int &read_gf,
   bool &finalize_topo, const std::string &xml_prefix)
{
   using namespace vtk_xml;

   static const char *erstr = "XML parsing error";

   // Create buffer beginning with xml_prefix, then read the rest of the stream
   std::vector<char> buf(xml_prefix.begin(), xml_prefix.end());
   std::istreambuf_iterator<char> eos;
   buf.insert(buf.end(), std::istreambuf_iterator<char>(input), eos);
   buf.push_back('\0'); // null-terminate buffer

   XMLDocument xml;
   xml.Parse(buf.data(), buf.size());
   if (xml.ErrorID() != XML_SUCCESS)
   {
      MFEM_ABORT("Error parsing XML VTK file.\n" << xml.ErrorStr());
   }

   const XMLElement *vtkfile = xml.FirstChildElement();
   MFEM_VERIFY(vtkfile, erstr);
   MFEM_VERIFY(StringCompare(vtkfile->Name(), "VTKFile"), erstr);
   const XMLElement *vtu = vtkfile->FirstChildElement();
   MFEM_VERIFY(vtu, erstr);
   MFEM_VERIFY(StringCompare(vtu->Name(), "UnstructuredGrid"), erstr);

   XMLDataReader data_reader(vtkfile, vtu);

   // Count the number of points and cells
   const XMLElement *piece = vtu->FirstChildElement();
   MFEM_VERIFY(StringCompare(piece->Name(), "Piece"), erstr);
   MFEM_VERIFY(piece->NextSiblingElement() == NULL,
               "XML VTK meshes with more than one Piece are not supported");
   int npts = piece->IntAttribute("NumberOfPoints");
   int ncells = piece->IntAttribute("NumberOfCells");

   // Read the points
   Vector points(3*npts);
   const XMLElement *pts_xml;
   for (pts_xml = piece->FirstChildElement();
        pts_xml != NULL;
        pts_xml = pts_xml->NextSiblingElement())
   {
      if (StringCompare(pts_xml->Name(), "Points"))
      {
         const XMLElement *pts_data = pts_xml->FirstChildElement();
         MFEM_VERIFY(pts_data->IntAttribute("NumberOfComponents") == 3,
                     "XML VTK Points DataArray must have 3 components");
         data_reader.Read(pts_data, points.GetData(), points.Size());
         break;
      }
   }
   if (pts_xml == NULL) { MFEM_ABORT(erstr); }

   // Read the cells
   Array<int> cell_data, cell_offsets(ncells), cell_types(ncells);
   const XMLElement *cells_xml;
   for (cells_xml = piece->FirstChildElement();
        cells_xml != NULL;
        cells_xml = cells_xml->NextSiblingElement())
   {
      if (StringCompare(cells_xml->Name(), "Cells"))
      {
         const XMLElement *cell_data_xml = NULL;
         for (const XMLElement *data_xml = cells_xml->FirstChildElement();
              data_xml != NULL;
              data_xml = data_xml->NextSiblingElement())
         {
            const char *data_name = data_xml->Attribute("Name");
            if (StringCompare(data_name, "offsets"))
            {
               data_reader.Read(data_xml, cell_offsets.GetData(), ncells);
            }
            else if (StringCompare(data_name, "types"))
            {
               data_reader.Read(data_xml, cell_types.GetData(), ncells);
            }
            else if (StringCompare(data_name, "connectivity"))
            {
               // Have to read the connectivity after the offsets, because we
               // don't know how many points to read until we have the offsets
               // (size of connectivity array is equal to the last offset), so
               // store the XML element pointer and read this data later.
               cell_data_xml = data_xml;
            }
         }
         MFEM_VERIFY(cell_data_xml != NULL, erstr);
         int cell_data_size = cell_offsets.Last();
         cell_data.SetSize(cell_data_size);
         data_reader.Read(cell_data_xml, cell_data.GetData(), cell_data_size);
         break;
      }
   }
   if (cells_xml == NULL) { MFEM_ABORT(erstr); }

   // Read the element attributes, which are stored as CellData named either
   // "material" or "attribute". We prioritize "material" over "attribute" for
   // backwards compatibility.
   Array<int> cell_attributes;
   bool found_attributes = false;
   for (const XMLElement *cell_data_xml = piece->FirstChildElement();
        cell_data_xml != NULL;
        cell_data_xml = cell_data_xml->NextSiblingElement())
   {
      const bool is_cell_data =
         StringCompare(cell_data_xml->Name(), "CellData");
      const bool is_material =
         StringCompare(cell_data_xml->Attribute("Scalars"), "material");
      const bool is_attribute =
         StringCompare(cell_data_xml->Attribute("Scalars"), "attribute");
      if (is_cell_data && (is_material || (is_attribute && !found_attributes)))
      {
         found_attributes = true;
         const XMLElement *data_xml = cell_data_xml->FirstChildElement();
         if (data_xml != NULL && StringCompare(data_xml->Name(), "DataArray"))
         {
            cell_attributes.SetSize(ncells);
            data_reader.Read(data_xml, cell_attributes.GetData(), ncells);
         }
      }

      // Simvascular specific attributes
      // read in modelfaceID for surface vtu file (sv meshing)
      if (StringCompare(cell_data_xml->Name(), "CellData")
          && StringCompare(cell_data_xml->Attribute("Scalars"), "ModelFaceID"))
      {
         const XMLElement *data_xml2;
         for (data_xml2 = cell_data_xml->FirstChildElement();
         data_xml2 != NULL;
         data_xml2 = data_xml2->NextSiblingElement())
         {
            if (data_xml2 != NULL && StringCompare(data_xml2->Attribute("Name"), "ModelFaceID"))
            {
               cell_attributes.SetSize(ncells);
               data_reader.Read(data_xml2, cell_attributes.GetData(), ncells);
               break;
            }
         }
      }

      // update 4-13-24
      if (StringCompare(cell_data_xml->Name(), "CellData")
          && StringCompare(cell_data_xml->Attribute("Normals"), "Normals"))
      {
         const XMLElement *data_xml2;
         for (data_xml2 = cell_data_xml->FirstChildElement();
         data_xml2 != NULL;
         data_xml2 = data_xml2->NextSiblingElement())
         {
            if (data_xml2 != NULL && StringCompare(data_xml2->Attribute("Name"), "ModelFaceID"))
            {
               cell_attributes.SetSize(ncells);
               data_reader.Read(data_xml2, cell_attributes.GetData(), ncells);
               break;
            }
         }
      }
   }

   CreateVTKMesh(points, cell_data, cell_offsets, cell_types, cell_attributes,
                 curved, read_gf, finalize_topo);
}
#endif