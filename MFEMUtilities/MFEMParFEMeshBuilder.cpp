#include "MFEMParFEMeshBuilder.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>

#include "../MFEMCore/MFEMMesh.h"
#include "../MFEMCore/MFEMLoad.h"
#include "../MFEMCore/MFEMMaterial.h"
#include "../MFEMCore/MFEMInitialCondition.h"
#include "../MFEMCore/MFEMBoundaryCondition.h"
#include "../MFEMImSolid/MFEMImmersedElasticDomain.h"
#include "mfem.hpp"

#ifdef MFEM_USE_MPI
MFEMParFEMeshBuilder::MFEMParFEMeshBuilder(FEModel* fem, MFEMMesh* mesh_) : m_fem(fem), m_serial_mesh(mesh_)
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

   // copy over surfaces
   int NSURF = mesh.Surfaces();
   serial_surf_list.SetSize(NSURF);
   for (int i=0; i<NSURF; ++i)
   {
      FESurface& surf = mesh.Surface(i);
      FEFacetSet *fs = surf.GetFacetSet();
      int surf_attr = m_serial_mesh->GetSurfaceAttribute(fs); assert(surf_attr > 0);
      
      serial_surf_list[surf_attr-1] = m_serial_femesh->CreateSurface(*fs);
      serial_surf_list[surf_attr-1]->Init();
   }
}

MFEMParFEMeshBuilder::~MFEMParFEMeshBuilder()
{
   if (m_serial_femesh) delete m_serial_femesh;
   if (m_febio_surf) delete m_febio_surf;
   ClearSerialSurfaces();
}

void MFEMParFEMeshBuilder::ClearSerialSurfaces()
{
   for (int i=0; i<serial_surf_list.Size(); ++i) 
   {
      if(serial_surf_list[i])
      {
         delete serial_surf_list[i];
         serial_surf_list[i] = nullptr;
      }
   }
}

void MFEMParFEMeshBuilder::BuildLocalFEBioMesh(mfem::ParFiniteElementSpace* dfes)
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
   global_node_ids_gf = -1.0;
   
   mfem::BaseKDTreeNodalProjection* map = new mfem::KDTreeNodalProjection<3>(global_node_ids_gf);
   map->Project(fe_node_ids_gf, 1e-12);
   delete map;
   
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
      node = global_node;
      node.SetID(i); // reset Node id, ensures local mfem id = febio id
	}

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

void MFEMParFEMeshBuilder::ParCopyDomain(FEDomain &dom, FEDomain &par_dom, mfem::ParFiniteElementSpace &pfes)
{
   auto sourceMesh = pfes.GetParMesh();
   par_dom.SetName(dom.GetName());

   // --- Elements --- //

   // Use a KDTree to identify the FEBio global element ids for each local mfem element
   // global_elem_ids_gf contains the febio global elem ids for each local mfem elem
   // after map->Project() below
   // Note this is done because after partitioning, the mfem global numbering of the elements
   // can change
   mfem::L2_FECollection l2_coll(0,3);
   mfem::FiniteElementSpace fes(m_submesh, &l2_coll, 1);
   mfem::Array<int> submesh_to_parent_elem_map = m_submesh->GetParentElementIDMap();
   mfem::GridFunction fe_elem_ids_gf(&fes);
   for (int i=0; i<m_submesh->GetNE(); ++i)
   {
      int parent_elem_id = submesh_to_parent_elem_map[i];
      fe_elem_ids_gf(i) = m_serial_mesh->GetFEBioElementID(parent_elem_id);
   }

   mfem::FiniteElementSpace mfes(sourceMesh, &l2_coll, 1);
   mfem::GridFunction global_elem_ids_gf(&mfes);
   global_elem_ids_gf = -1.0;
   
   mfem::BaseKDTreeNodalProjection* map = new mfem::KDTreeNodalProjection<3>(global_elem_ids_gf);
   map->Project(fe_elem_ids_gf, 1e-12);
   delete map;

   int NE = pfes.GetNE();
   par_dom.Create(NE,dynamic_cast<FESolidDomain*>(&dom)->GetElementSpec());
   mfem::Array<int> dofs;

   // loop over local elements in sourceMesh
   std::map<int,int> febio_to_local_elem_map; // map local mfem elem id to global febio elem id
   for (int i=0; i<NE; ++i)
   {
      int glob_febio_elem_id = global_elem_ids_gf(i);
      MFEM_ASSERT(glob_febio_elem_id >= 0, "Error: element not found");

      // determine which domain this element belongs to
      MFEMImmersedElasticDomain *immersed_dom = nullptr;
      FESolidElement* el = nullptr;
      for (int j=0; j<m_serial_femesh->Domains(); ++j)
      {
         FEDomain& test_dom = m_serial_femesh->Domain(j);
         immersed_dom = dynamic_cast<MFEMImmersedElasticDomain*>(&test_dom);
         if (immersed_dom == nullptr) { continue; }
         el = dynamic_cast<FESolidElement*>(immersed_dom->FindElementFromID(glob_febio_elem_id));
         if (el != nullptr)
         {
            // found the element in the immersed domain
            break;
         }
      }

      if (el == nullptr)
      {
         MFEM_ABORT("Error: Could not find element in immersed domain, elem number = " + std::to_string(glob_febio_elem_id));
      }

      // copy the element data to new domain element
      FESolidElement* pel = &dynamic_cast<FESolidDomain*>(&par_dom)->Element(i);
      *pel = *el;
      pel->SetID(i+1); // set local febio element id = mfem element id + 1

      // reset some data in pel
      // reset the material point pointer in each Gauss point
      for (int j=0; j<pel->GaussPoints(); ++j)
      {
         FEMaterialPoint* pmp = pel->GetMaterialPoint(j);
         if (pmp) pmp->m_elem = pel; // reset the element pointer in material point
      }

      febio_to_local_elem_map[glob_febio_elem_id] = i;

      // update element node list
      pfes.GetElementDofs(i,dofs);
      
      MFEM_ASSERT(pel->Nodes() == dofs.Size(), "Error: Number of nodes in FEBio element does not match number of dofs in mfem element");
      for (int j=0; j<pel->Nodes(); ++j)
      {
         pel->m_node[j] = dofs[j];
      }     
      // set the mesh partition
      pel->SetMeshPartition(&par_dom);
   }

   // set the mesh of this domain to the local fe mesh
   par_dom.SetMesh(&m_fem->GetMesh());

   // Surface Elements //

   // Use a KDTree to identify the global surface element ids for each local mfem surface element
   // global_surf_elem_ids_gf contains the global elem ids for each local mfem surface elem
   // after map->Project() below

   // Extract global bdr from the mfem serial mesh
   auto global_bdr_submesh = new mfem::SubMesh(mfem::SubMesh::CreateFromBoundary(*m_serial_mesh,m_serial_mesh->bdr_attributes));   
   // maps global bdrsubmesh element to global surface element
   mfem::Array<int> bdr_to_global_map = global_bdr_submesh->GetParentElementIDMap();

   // Extract solid bdr from global bdr submesh
   // Note this will ignore extra bdr attributes generated by mfem submesh when extracting the solid boundary
   if (global_bdr_submesh->bdr_attributes.Size() == 0)
   {
      // append a single bdr attribute if none exist
      // this ensures CreateFromDomain() below works correctly
      // its not used for anything
      global_bdr_submesh->bdr_attributes.Append(1);
   }
   auto solid_bdr_submesh = new mfem::SubMesh(mfem::SubMesh::CreateFromDomain(*global_bdr_submesh,m_submesh->bdr_attributes));
   
   // no solid boundary elements found (no bcs applied to solid boundary),
   // then we can skip creating the surface elements
   if (solid_bdr_submesh->GetNE() == 0)
   {      
      delete global_bdr_submesh;
      delete solid_bdr_submesh;
   }
   else
   {
      mfem::Array<int> solid_bdr_to_bdr_map = solid_bdr_submesh->GetParentElementIDMap();
      mfem::L2_FECollection l2_coll_bdr(0,2);
      mfem::FiniteElementSpace solid_bdr_fes(solid_bdr_submesh, &l2_coll_bdr, 1);
      
      // Loop over the solid boundary elements in the submesh
      mfem::GridFunction fe_surf_elem_ids_gf(&solid_bdr_fes);
      for (int i=0; i<solid_bdr_submesh->GetNE(); ++i)
      {
         int bdr_id = solid_bdr_to_bdr_map[i];
         int global_bdr_id = bdr_to_global_map[bdr_id];

         fe_surf_elem_ids_gf(i) = global_bdr_id;
      }
      
      // create a bdr_submesh of the sourceMesh (parallel partitioned mfem mesh)
      auto bdr_sourceMesh = new mfem::SubMesh(mfem::SubMesh::CreateFromBoundary(*sourceMesh,sourceMesh->bdr_attributes));
      mfem::FiniteElementSpace bdr_mfes(bdr_sourceMesh, &l2_coll_bdr, 1);
      
      mfem::GridFunction global_surf_elem_ids_gf(&bdr_mfes);
      global_surf_elem_ids_gf = -1.0; // initialize to -1 to indicate bdr elements in immersed solid that do not have a BC assigned

      mfem::BaseKDTreeNodalProjection* map_surf = new mfem::KDTreeNodalProjection<3>(global_surf_elem_ids_gf);
      map_surf->Project(fe_surf_elem_ids_gf, 1e-12);

      // free memory
      delete map_surf;
      delete bdr_sourceMesh;
      delete global_bdr_submesh;
      delete solid_bdr_submesh;

      // Loop over local surface elements and create the local FEBio surface elements
      m_febio_surf = new FESurface(m_fem);
      m_febio_surf->Create(pfes.GetNBE());
      mfem::Array<int> vdofs;
      for (int i=0; i<pfes.GetNBE(); ++i)
      {
         // get the mfem global surface element id
         int global_mfem_surf_elem_id = global_surf_elem_ids_gf(i);
         if (global_mfem_surf_elem_id < 0)
         {
            // this surface element does not have a boundary condition applied
            // so we skip it

            FESurfaceElement* surf_el = &m_febio_surf->Element(i);
            surf_el->SetID(-1); // mark surface element as inactive

            continue;
         }

         // get the febio surface element
         int bdr_attr = sourceMesh->GetBdrAttribute(i);
         int febio_surf_id = m_serial_mesh->GetFEBioSurfaceElementID(bdr_attr,global_mfem_surf_elem_id);
         FESurface& surf = *serial_surf_list[bdr_attr - 1];
         MFEM_ASSERT(surf.IsActive(), "Error: the current surface is not active");
         FESurfaceElement* el = &surf.Element(febio_surf_id);
         MFEM_ASSERT(el != nullptr, "Error: Could not find surface element in febio mesh, surface element = " + std::to_string(global_mfem_surf_elem_id));
         
         // copy data from febio surface element to new local surface element
         FESurfaceElement* pel = &m_febio_surf->Element(i);
         *pel = *el;
         pel->SetID(i+1); // set local febio surface element id = mfem surface element id + 1
         
         // reset dofs to match mfem
         pfes.GetBdrElementDofs(i, vdofs);
         MFEM_ASSERT(vdofs.Size() == pel->Nodes(), "Error: Number of dofs in febio surface element does not match number of dofs in mfem surface element");
         for (int j=0; j<vdofs.Size(); ++j)
         {
            pel->m_node[j] = vdofs[j];
         }

         // create material point data
         pel->ClearData();
         int nint = pel->GaussPoints();
         for (int n = 0; n < nint; ++n)
         {
            FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(new FESurfaceMaterialPoint);
            assert(pt);
            pel->SetMaterialPointData(pt, n);
         }
      }

      vec3d re[FEElement::MAX_NODES];
      FEMesh &mesh = m_fem->GetMesh();
      // initialize material points of surface elements
      for (int i = 0; i < m_febio_surf->Elements(); ++i)
      {
         FESurfaceElement& el = m_febio_surf->Element(i);
         if (el.GetID() < 0) continue; // skip inactive surface elements

         // m_febio_surf->NodalCoordinates(el, re);
         for (int ii = 0; ii < el.Nodes(); ++ii) re[ii] = mesh.Node(el.m_lnode[ii]).m_rt;

         int nint = el.GaussPoints();
         int neln = el.Nodes();
         for (int n = 0; n < nint; ++n)
         {
            FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(el.GetMaterialPoint(n));
            if (pt == nullptr) MFEM_ABORT("Error: Could not create surface material point");

            // initialize some material point data
            double* H = el.H(n);
            vec3d rn(0, 0, 0);
            for (int j = 0; j < neln; ++j)
            {
               rn += re[j]*H[j];
            }

            pt->m_r0 = rn;

            // calculate initial surface tangents
            double* Gr = el.Gr(n);
            double* Gs = el.Gs(n);

            vec3d dxr(0, 0, 0), dxs(0, 0, 0);
            for (int i = 0; i < neln; ++i)
            {
               dxr += re[i] * Gr[i];
               dxs += re[i] * Gs[i];
            }
            pt->dxr = dxr;
            pt->dxs = dxs;

            // initialize the other material point data
            pt->Init();
         }
      }

      // Add the combined surface to the FEMesh
      m_fem->GetMesh().AddSurface(m_febio_surf);
   }
}

#endif