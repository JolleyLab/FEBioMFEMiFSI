#include "ParFEModel.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>

#include "MFEMMesh.h"
#include "MFEMLoad.h"
#include "MFEMInitialCondition.h"
#include "MFEMBoundaryCondition.h"
#include "MFEMImmersedElasticDomain.h"
#include "mfem.hpp"

#ifdef MFEM_USE_MPI
ParFEModel::ParFEModel(FEModel* fem) : m_fem(fem)
{
   m_serial_femesh = new FEMesh(fem);
   m_serial_femesh->CopyFrom(fem->GetMesh());

   // loop over domains and set materials
   FEMesh &mesh = m_fem->GetMesh();
   int NDOM = mesh.Domains();
   for (int i=0; i<NDOM; ++i)
   {
      FEDomain& dom = mesh.Domain(i);
      FEDomain& target_dom = m_serial_femesh->Domain(i);
      FEMaterial* pmat = dom.GetMaterial();
      if (pmat) { target_dom.SetMaterial(pmat); }
      else assert(false);
   }
}

ParFEModel::~ParFEModel()
{
   if (m_serial_femesh) delete m_serial_femesh;
}

void ParFEModel::SetMFEMSerialMesh(MFEMMesh* mesh)
{
   m_serial_mesh = mesh;
}

void ParFEModel::SetupBCs()
{
   // --- boundary conditions ---

	int NDC = m_fem->BoundaryConditions();
	for (int i=0; i<NDC; ++i)
	{
		FEBoundaryCondition* pbc = m_fem->BoundaryCondition(i);
      MFEMBoundaryCondition* m_bc = dynamic_cast<MFEMBoundaryCondition*>(pbc);
      assert(m_bc);

      FESurface &surf = *m_bc->GetSurface();
      FEFacetSet* set = surf.GetFacetSet();
      int id = m_serial_mesh->GetSurfaceAttribute(set); assert(id > 0);
      m_bc->SetBoundaryID(id);
	}

   // --- model loads ---

   int NML = m_fem->ModelLoads();
	for (int i=0; i<NML; ++i)
	{
		FEModelLoad* pml = m_fem->ModelLoad(i);

      // nodal loads
      MFEMNodalLoad* nl = dynamic_cast<MFEMNodalLoad*>(pml);
      // surface loads
		MFEMSurfaceLoad* sl = dynamic_cast<MFEMSurfaceLoad*>(pml);
      // body loads
      MFEMBodyLoad* bl = dynamic_cast<MFEMBodyLoad*>(pml);

      if (nl) 
      {
         MFEM_ABORT("TODO");
         // nl->CopyFrom(pml);
      }
      else if (sl)
      {
         FESurface &surf = sl->GetSurface();
         FEFacetSet* set = surf.GetFacetSet();
         int id = m_serial_mesh->GetSurfaceAttribute(set); assert(id > 0);
         sl->SetBoundaryID(id);
      } 
      else if (bl) 
      {
         std::vector<int> dom_list;
         for (int j=0; j<bl->Domains(); ++j)
         {
            FEDomain* dom = bl->Domain(j);
            int id = m_serial_mesh->GetDomainAttribute(dom); assert(id > 0);
            dom_list.emplace_back(id);
         }   
         bl->SetDomainIDList(dom_list);
      }
      else MFEM_ABORT("Error: unknown load type");
	}

   // // --- initial conditions ---
   
   // int NIC = m_fem->InitialConditions();
	// for (int i=0; i<NIC; ++i)
	// {
	// 	FEInitialCondition* pic = m_fem->InitialCondition(i);
   //    MFEMInitialCondition* mic = dynamic_cast<MFEMInitialCondition*>(pic);
   //    assert(mic);
	// 	mic->CopyFrom(*pic);
	// }
}

void ParFEModel::BuildLocalFEBioMesh(mfem::ParFiniteElementSpace* dfes)
{
   // --- Mesh ---

	// copy the mesh data
	FEMesh& mesh = m_fem->GetMesh(); // mesh to rebuild from mfem ParMesh
   auto sourceMesh = dfes->GetParMesh();
   mfem::GridFunction &Nodes = *sourceMesh->GetNodes();
   int NV = sourceMesh->GetNV();

   // --- Vertices --- //

   // Use a KDTree to identify the FEBio global node ids for each local mfem node
   // global_node_ids_gf contains the febio global node ids for each local mfem node
   // after map->Project() below
   mfem::H1_FECollection h1_coll(1,3);
   mfem::FiniteElementSpace fes(m_submesh, &h1_coll, 1);
   mfem::Array<int> array_submesh_to_parent_vertices = m_submesh->GetParentVertexIDMap();
   mfem::GridFunction fe_node_ids_gf(&fes);
   for (int i=0; i<m_submesh->GetNV(); ++i)
   {
      // need to map submesh vertex id back to parent mesh
      fe_node_ids_gf(i) = array_submesh_to_parent_vertices[i];
   }

   mfem::FiniteElementSpace mfes(sourceMesh, &h1_coll, 1);
   mfem::GridFunction global_node_ids_gf(&mfes);
   global_node_ids_gf = 0.0;
   
   mfem::BaseKDTreeNodalProjection* map = new mfem::KDTreeNodalProjection<3>(global_node_ids_gf);
   map->Project(fe_node_ids_gf, 1e-12);
   delete map;
   // for (int i=0; i<sourceMesh->GetNV(); ++i)
   // {
   //    mfem::out << "Node[" << i << "] = " << "FEBio Node[" << global_node_ids_gf(i)+1 << "]\n";
   // }   
   
	mesh.CreateNodes(NV);
   // loop over vertices in sourceMesh and update m_Node[] list in mesh
	for (int i=0; i<NV; ++i)
	{
      // get febio global node id
      int global_febio_id = global_node_ids_gf(i);

      // get the node from FEBio source mesh
      MFEM_ASSERT(global_febio_id < m_serial_femesh->Nodes(), "Error: global node index out of range");
      MFEM_ASSERT(global_febio_id >= 0, "Error: global node index out of range");
      FENode& global_node = m_serial_femesh->Node(global_febio_id);

      // copy data to new node
      FENode &node = mesh.Node(i);
      // int nid = node.GetID();
      node = global_node;
      node.SetID(i); // reset Node id, ensures local mfem id = febio id

      // //print out nodes
      // std::cout << "--------------------------------\n";
      // std::cout << "Local mfem id = " << i << ", Global FEBio id = " << global_febio_id + 1 << "\n";
      // std::cout << "MFEM Node = (" 
      //          << Nodes(i*3 + 0) << ", " 
      //          << Nodes(i*3 + 1) << ", " 
      //          << Nodes(i*3 + 2) << ")\n";
      // std::cout << "FEBio Node = (" 
      //          << global_node.m_r0.x << ", " 
      //          << global_node.m_r0.y << ", " 
      //          << global_node.m_r0.z << ")\n";

      // // check if any differ
      // double diffx = Nodes(i*3 + 0) - global_node.m_r0.x;
      // double diffy = Nodes(i*3 + 1) - global_node.m_r0.y;
      // double diffz = Nodes(i*3 + 2) - global_node.m_r0.z;
      // if (std::abs(diffx) > 1e-10 || 
      //     std::abs(diffy) > 1e-10 || 
      //     std::abs(diffz) > 1e-10)
      // {
      //    MFEM_ABORT("Node differs!");
      // }

	}

   // MFEM_ABORT("HERE");

   // --- Domains --- //

	// B. domains (only consider solid domains for now)
   // In building the mfem serial solid mesh, we have already
   // combined the immersed solid domains into a single domain
   // so we can just copy the mfem serial solid mesh to a single domain

   // TODO: Add support for immersed solids with different material types
   // would need to build each domain separately!
   
   // find the immersed domain
   int NDOM = m_serial_femesh->Domains();
   MFEMImmersedElasticDomain* mid = nullptr;
	for (int i=0; i<NDOM; ++i)
	{
		FEDomain& dom = m_serial_femesh->Domain(i);
      mid = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (mid) { break; }
   }
   const char* sz = mid->GetTypeStr();

   // create a new domain
   FEDomain* pd = nullptr;
   switch (mid->Class())
   {
   case FE_DOMAIN_SOLID   : pd = fecore_new<FESolidDomain   >(sz, m_fem); break;
   // case FE_DOMAIN_SHELL   : pd = fecore_new<FEShellDomain   >(sz, this); break;
   // case FE_DOMAIN_BEAM    : pd = fecore_new<FEBeamDomain    >(sz, this); break;
   // case FE_DOMAIN_2D      : pd = fecore_new<FEDomain2D      >(sz, this); break;
   // case FE_DOMAIN_DISCRETE: pd = fecore_new<FEDiscreteDomain>(sz, this); break;
   }
   assert(pd);
   pd->SetMaterial(mid->GetMaterial()); // here we assume one material type for the immersed solid domain

   // create the domain to only include elements contained in this partition
   ParCopyDomain(*mid,*pd,*dfes);

   // re-init the domain
   dynamic_cast<FESolidDomain*>(pd)->Init();   

   // add it to the mesh
   mesh.AddDomain(pd);

   // clear the copy of febio mesh
   m_serial_femesh->Clear();
}

void ParFEModel::ParCopyDomain(FEDomain &dom, FEDomain &par_dom, mfem::ParFiniteElementSpace &pfes)
{
   auto sourceMesh = pfes.GetParMesh();
   par_dom.SetName(dom.GetName());

   // --- Elements --- //

   // Use a KDTree to identify the FEBio global element ids for each local mfem element
   // global_elem_ids_gf contains the febio global elem ids for each local mfem elem
   // after map->Project() below
   mfem::L2_FECollection l2_coll(0,3);
   mfem::FiniteElementSpace fes(m_submesh, &l2_coll, 1);
   mfem::Array<int> submesh_to_parent_elem_map = m_submesh->GetParentElementIDMap();
   mfem::GridFunction fe_elem_ids_gf(&fes);
   for (int i=0; i<m_submesh->GetNE(); ++i)
   {
      fe_elem_ids_gf(i) = submesh_to_parent_elem_map[i] + 1; // febio element ids are 1-based, mfem is 0-based
   }

   mfem::FiniteElementSpace mfes(sourceMesh, &l2_coll, 1);
   mfem::GridFunction global_elem_ids_gf(&mfes);
   global_elem_ids_gf = 0.0;
   
   mfem::BaseKDTreeNodalProjection* map = new mfem::KDTreeNodalProjection<3>(global_elem_ids_gf);
   map->Project(fe_elem_ids_gf, 1e-12);
   delete map;

   int NE = pfes.GetNE();
   par_dom.Create(NE,dynamic_cast<FESolidDomain*>(&dom)->GetElementSpec());
   mfem::Array<int> dofs;

   // mfem::Array<int> submesh_to_parent_elem_map = m_submesh->GetParentElementIDMap();
   // mfem::Vector elem_center(3);

   // loop over local elements in sourceMesh
   for (int i=0; i<NE; ++i)
   {
      // // get global mfem submesh element id
      // int glob_submesh_elem_id = sourceMesh->GetGlobalElementNum(i);

      // // // get global febio element id
      // // // NOTE: febio is 1-based, mfem is 0-based for element numbering
      // // // global elem id (febio is 1-based, mfem is 0-based)
      // // // int glob_febio_elem_id = submesh_to_parent_elem_map[glob_submesh_elem_id] + 1;
      // int glob_febio_elem_id2 = m_serial_mesh->GetFEBioElementID(submesh_to_parent_elem_map[glob_submesh_elem_id]);

      // // sourceMesh->GetElementCenter(i,elem_center);
      // // std::vector<double> elem_center_std(3);
      // // elem_center_std[0] = elem_center(0);
      // // elem_center_std[1] = elem_center(1);
      // // elem_center_std[2] = elem_center(2);
      // // int glob_febio_elem_id = m_serial_mesh->GetFEBioElementID(elem_center_std);
      // // std::cout << "Local mfem element id = " << i << ", Global mfem element id = " << glob_submesh_elem_id << ", Global FEBio elem id = " << glob_febio_elem_id << "\n";
      // // std::cout << "Element center = (" 
      // //           << elem_center_std[0] << ", " 
      // //           << elem_center_std[1] << ", " 
      // //           << elem_center_std[2] << ")\n";

      int glob_febio_elem_id = global_elem_ids_gf(i);

      // determine which domain this element belongs to
      MFEMImmersedElasticDomain *immersed_dom = nullptr;
      FESolidElement* el = nullptr;
      for (int j=0; j<m_serial_femesh->Domains(); ++j)
      {
         FEDomain& test_dom = m_serial_femesh->Domain(j);
         immersed_dom = dynamic_cast<MFEMImmersedElasticDomain*>(&test_dom);
         if (immersed_dom == nullptr) { continue; }
         // FESolidElement& el = immersed_dom->Element(0);
         // int elem_lower = el.GetID();
         // el = immersed_dom->Element(immersed_dom->Elements()-1);
         // int elem_upper = el.GetID();
         
         // if (glob_febio_elem_id >= elem_lower && glob_febio_elem_id <= elem_upper)
         // {
         //    break;
         // }

         el = dynamic_cast<FESolidElement*>(immersed_dom->FindElementFromID(glob_febio_elem_id));
         if (el != nullptr)
         {
            // found the element in the immersed domain
            break;
         }
      }
      MFEM_ASSERT(el != nullptr, "Error: Could not find element in immersed domain, elem number = " + std::to_string(glob_febio_elem_id));

      // copy the element data to new domain element
      FESolidElement* pel = &dynamic_cast<FESolidDomain*>(&par_dom)->Element(i);
      *pel = *el;
      pel->SetID(i); // set local febio element id = mfem element id

      // update element node list
      pfes.GetElementDofs(i,dofs);
      
      
      // std::cout << "febio num nodes = " << pel->Nodes() << ", mfem num dofs = " << dofs.Size() << std::endl;
     
      // auto fe = pfes.GetFE(i);
      // auto gtype = fe->GetGeomType();
      // std::string mfem_type;
      // switch (gtype)
      // {
      //    case mfem::Geometry::POINT: mfem_type = "POINT"; break;
      //    case mfem::Geometry::SEGMENT: mfem_type = "SEGMENT"; break;
      //    case mfem::Geometry::TRIANGLE: mfem_type = "TRIANGLE"; break;
      //    case mfem::Geometry::SQUARE: mfem_type = "QUADRILATERAL"; break;
      //    case mfem::Geometry::TETRAHEDRON: mfem_type = "TETRAHEDRON"; break;
      //    case mfem::Geometry::CUBE: mfem_type = "HEXAHEDRON"; break;
      //    case mfem::Geometry::PRISM: mfem_type = "PRISM"; break;
      //    case mfem::Geometry::PYRAMID: mfem_type = "PYRAMID"; break;
      //    default: mfem_type = "UNKNOWN"; break;
      // }

      // std::string febio_type;
      // switch (pel->Type())
      // {
      //    // Hexahedral elements
      //    case FE_HEX8G8:
      //    case FE_HEX8RI:
      //    case FE_HEX8G1:
      //    case FE_HEX20G8:
      //    case FE_HEX20G27:
      //    case FE_HEX27G27:
      //       febio_type = "HEXAHEDRON"; 
      //       break;
            
      //    // Tetrahedral elements
      //    case FE_TET4G1:
      //    case FE_TET4G4:
      //    case FE_TET5G4:
      //    case FE_TET10G1:
      //    case FE_TET10G4:
      //    case FE_TET10G8:
      //    case FE_TET10GL11:
      //    case FE_TET10G4RI1:
      //    case FE_TET10G8RI4:
      //    case FE_TET15G4:
      //    case FE_TET15G8:
      //    case FE_TET15G11:
      //    case FE_TET15G15:
      //    case FE_TET15G15RI4:
      //    case FE_TET20G15:
      //       febio_type = "TETRAHEDRON"; 
      //       break;
            
      //    // Pentahedral/Wedge elements
      //    case FE_PENTA6G6:
      //    case FE_PENTA15G8:
      //    case FE_PENTA15G21:
      //       febio_type = "PRISM"; 
      //       break;
            
      //    // Pyramid elements
      //    case FE_PYRA5G8:
      //    case FE_PYRA13G8:
      //       febio_type = "PYRAMID"; 
      //       break;
            
      //    // Shell elements
      //    case FE_SHELL_QUAD4G4:
      //    case FE_SHELL_QUAD4G8:
      //    case FE_SHELL_QUAD4G12:
      //    case FE_SHELL_QUAD8G18:
      //    case FE_SHELL_QUAD8G27:
      //    case FE_SHELL_TRI3G3:
      //    case FE_SHELL_TRI3G6:
      //    case FE_SHELL_TRI3G9:
      //    case FE_SHELL_TRI6G14:
      //    case FE_SHELL_TRI6G21:
      //       febio_type = "SHELL"; 
      //       break;
            
      //    // Beam elements
      //    case FE_BEAM2G1:
      //    case FE_BEAM2G2:
      //    case FE_BEAM3G2:
      //       febio_type = "BEAM"; 
      //       break;
            
      //    // Discrete elements
      //    case FE_DISCRETE:
      //       febio_type = "DISCRETE"; 
      //       break;
            
      //    // Truss elements
      //    case FE_TRUSS:
      //       febio_type = "TRUSS"; 
      //       break;
            
      //    default: 
      //       febio_type = "UNKNOWN"; 
      //       break;
      // }

      // std::cout << "febio element type = " << febio_type << ", mfem element type = " << mfem_type << std::endl;
      // std::cout << "Local mfem element id = " << i << ", Global mfem element id = " << glob_submesh_elem_id << ", Global FEBio elem id = " << glob_febio_elem_id << "\n";
      // // print out mfem node coordinates
      // // print out nodal coordinates of each node on the mfem mesh
      // auto Nodes1 = pfes.GetFE(i)->GetNodes();
      // auto eltrans1 = pfes.GetElementTransformation(i);
      // std::cout << "MFEM Nodes = " << "\n";
      // mfem::Vector curr_pos1(3);
      // for(int j=0; j<Nodes1.GetNPoints(); ++j)
      // {
      //    auto ip_node = Nodes1.IntPoint(j);
      //    eltrans1->SetIntPoint(&ip_node);
      //    eltrans1->Transform(ip_node, curr_pos1);
      //    std::cout << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos1(0) << "," << curr_pos1(1) << "," << curr_pos1(2) << ")" << std::endl;
      // }    
      
      MFEM_ASSERT(pel->Nodes() == dofs.Size(), "Error: Number of nodes in FEBio element does not match number of dofs in mfem element");
      for (int j=0; j<pel->Nodes(); ++j)
      {
         pel->m_node[j] = dofs[j];
      }     

      // // print out febio node coordinates
      // std::cout << "\n";
      // std::cout << "--------------------------------\n";
      // std::cout << "Local mfem element id = " << i << ", Global mfem element id = " << glob_submesh_elem_id << ", Global FEBio elem id = " << glob_febio_elem_id << "\n";
      // // std::cout << "Element id = " << glob_febio_elem_id << std::endl;
      // std::cout << "FEBio Nodes = " << "\n";
      // for (int j=0; j<pel->Nodes(); ++j)
      // {
      //    FENode& node = m_fem->GetMesh().Node(pel->m_node[j]);
      //    std::cout << std::setprecision(8) << "Node[" << j << "]=(" << node.m_r0(0) << "," << node.m_r0(1) << "," << node.m_r0(2) << ")" << std::endl;
      // }

      // // print out mfem node coordinates
      // // print out nodal coordinates of each node on the mfem mesh
      // auto Nodes = pfes.GetFE(i)->GetNodes();
      // auto eltrans = pfes.GetElementTransformation(i);
      // std::cout << "MFEM Nodes = " << "\n";
      // mfem::Vector curr_pos(3);
      // for(int j=0; j<Nodes.GetNPoints(); ++j)
      // {
      //    auto ip_node = Nodes.IntPoint(j);
      //    eltrans->SetIntPoint(&ip_node);
      //    eltrans->Transform(ip_node, curr_pos);
      //    std::cout << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
      // }     
      
      // // check if any node differs
      // for(int j=0; j<pel->Nodes(); ++j)
      // {
      //    FENode& node = m_fem->GetMesh().Node(pel->m_node[j]);
      //    auto ip_node = Nodes.IntPoint(j);
      //    eltrans->SetIntPoint(&ip_node);
      //    eltrans->Transform(ip_node, curr_pos);

      //    double diffx = node.m_r0.x - curr_pos(0);
      //    double diffy = node.m_r0.y - curr_pos(1);
      //    double diffz = node.m_r0.z - curr_pos(2);
      //    if (fabs(diffx) > 1e-10 || fabs(diffy) > 1e-10 || fabs(diffz) > 1e-10)
      //    {
      //       std::cout << "Node[" << j << "] differs: FEBio=(" << node.m_r0.x << "," << node.m_r0.y << "," << node.m_r0.z << "), "
      //                 << "MFEM=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
      //       MFEM_ABORT("Node coordinates differ between FEBio and MFEM meshes!");
      //    }
      // }

      // set the mesh partition
      pel->SetMeshPartition(&par_dom);
   }

   // MFEM_ABORT("HERE");

   // set the mesh of this domain to the local fe mesh
   par_dom.SetMesh(&m_fem->GetMesh());
}

#endif