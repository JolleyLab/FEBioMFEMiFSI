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
   mfem::FiniteElementSpace fes(m_serial_mesh, &h1_coll, 1);
   mfem::GridFunction fe_node_ids_gf(&fes);
   for (int i=0; i<m_serial_mesh->GetNV(); ++i)
   {
      fe_node_ids_gf(i) = i;
   }

   mfem::FiniteElementSpace mfes(sourceMesh, &h1_coll, 1);
   mfem::GridFunction global_node_ids_gf(&mfes);
   global_node_ids_gf = 0.0;
   
   mfem::BaseKDTreeNodalProjection* map = new mfem::KDTreeNodalProjection<3>(global_node_ids_gf);
   map->Project(fe_node_ids_gf, 1e-10);
   delete map;
   // for (int i=0; i<sourceMesh->GetNV(); ++i)
   // {
   //    mfem::out << "Node[" << i << "] = " << "FEBio Node[" << global_node_ids_gf(i)+1 << "]\n";
   // }   
   
	mesh.CreateNodes(NV);
   // loop over vertices in sourceMesh and update m_Node[] list
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
      // mfem::out << "--------------------------------\n";
      // mfem::out << "Local mfem id = " << i << ", Global FEBio id = " << global_febio_id + 1 << "\n";
      // mfem::out << "MFEM Node = (" 
      //          << Nodes(i*3 + 0) << ", " 
      //          << Nodes(i*3 + 1) << ", " 
      //          << Nodes(i*3 + 2) << ")\n";
      // mfem::out << "FEBio Node = (" 
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

   int NE = pfes.GetNE();
   par_dom.Create(NE,dynamic_cast<FESolidDomain*>(&dom)->GetElementSpec());
   mfem::Array<int> dofs;
   mfem::Array<int> submesh_to_parent_elem_map = m_submesh->GetParentElementIDMap();
   // loop over local elements in sourceMesh
   for (int i=0; i<NE; ++i)
   {
      // get global mfem submesh element id
      int glob_submesh_elem_id = sourceMesh->GetGlobalElementNum(i);

      // get global febio element id
      // NOTE: febio is 1-based, mfem is 0-based for element numbering
      // global elem id (febio is 1-based, mfem is 0-based)
      int glob_febio_elem_id = submesh_to_parent_elem_map[glob_submesh_elem_id] + 1;

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

      // copy the element data to new domain element
      FESolidElement* pel = &dynamic_cast<FESolidDomain*>(&par_dom)->Element(i);
      *pel = *el;
      pel->SetID(i); // set local febio element id = mfem element id

      // update element node list
      pfes.GetElementDofs(i,dofs);
      for (int j=0; j<pel->Nodes(); ++j)
      {
         pel->m_node[j] = dofs[j];
      }     
      
      // // print out febio node coordinates
      // mfem::out << "\n";
      // mfem::out << "--------------------------------\n";
      // mfem::out << "Local mfem element id = " << i << ", Global mfem element id = " << glob_submesh_elem_id << ", Global FEBio elem id = " << glob_febio_elem_id << "\n";
      // // mfem::out << "Element id = " << glob_febio_elem_id << std::endl;
      // mfem::out << "FEBio Nodes = " << "\n";
      // for (int j=0; j<pel->Nodes(); ++j)
      // {
      //    FENode& node = m_fem->GetMesh().Node(pel->m_node[j]);
      //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << node.m_r0(0) << "," << node.m_r0(1) << "," << node.m_r0(2) << ")" << std::endl;
      // }

      // // print out mfem node coordinates
      // // print out nodal coordinates of each node on the mfem mesh
      // auto Nodes = pfes.GetFE(i)->GetNodes();
      // auto eltrans = pfes.GetElementTransformation(i);
      // mfem::out << "MFEM Nodes = " << "\n";
      // mfem::Vector curr_pos(3);
      // for(int j=0; j<Nodes.GetNPoints(); ++j)
      // {
      //    auto ip_node = Nodes.IntPoint(j);
      //    eltrans->SetIntPoint(&ip_node);
      //    eltrans->Transform(ip_node, curr_pos);
      //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
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
      //       mfem::out << "Node[" << j << "] differs: FEBio=(" << node.m_r0.x << "," << node.m_r0.y << "," << node.m_r0.z << "), "
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