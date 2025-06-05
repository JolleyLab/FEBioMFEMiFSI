#include "MFEMMesh.h"
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESurface.h>
#include "mfem.hpp"
#include <math.h>

MFEMMesh::MFEMMesh(int nodes, int elems, int bndrelems) : mfem::Mesh(3, nodes, elems, bndrelems)
{

}

MFEMMesh* CreateMFEMMesh3D(FEMesh& febioMesh)
{
	// We count the nr of boundary elements, by adding up all faces from the facet sets. 
	// Note that this may not count all the elements on the mesh boundary, but it doesn't look like that is a problem.
	// This also assumes that none of the facet sets overlap.
	int nbndrElems = 0;
	for (int n = 0; n < febioMesh.FacetSets(); ++n)
	{
		FEFacetSet& surf = febioMesh.FacetSet(n);
		int faces = surf.Faces();
		nbndrElems += faces;
	}

	int nodes = febioMesh.Nodes();
	int elems = febioMesh.Elements();
	MFEMMesh* newMesh = new MFEMMesh(nodes, elems, nbndrElems);
	for (int i = 0; i < nodes; ++i)
	{
		FENode& node = febioMesh.Node(i);
		vec3d& pos = node.m_r0;
		double coords[3] = { pos.x, pos.y, pos.z };
		newMesh->AddVertex(coords);
	}

	int indices[FEElement::MAX_NODES] = { 0 };
   int elem_counter = 0;
   mfem::Vector elem_center(3);
   std::vector<double> elem_center_std(3);
   int new_elem_num = 0;
	for (int n = 0; n < febioMesh.Domains(); ++n)
	{
		FEDomain& dom = febioMesh.Domain(n);
		int elemsInDomain = dom.Elements();
      newMesh->AssignDomainAttribute(&dom, n + 1);
		for (int i = 0; i < elemsInDomain; i++) {
         std::fill(indices, indices + FEElement::MAX_NODES, 0); // reset indices to all 0
			FEElement& el = dom.ElementRef(i);
			for (int j = 0; j < el.Nodes(); j++) {
				indices[j] = el.m_node[j];
			}

         // // Set the FEBio element id for the mfem element
         // newMesh->SetFEBioElementID(elem_counter, el.GetID());
         // // if (elem_counter == 8482 + 620 || elem_counter == 8482 + 104)
         // // {
         // //    int val = 0;
         // // } 
         // elem_counter++;

			switch (el.Shape())
			{
         case ET_TET4: 
         case ET_TET10:
            new_elem_num = newMesh->AddTet(indices, n + 1); break;
         case ET_HEX8: 
         case ET_HEX27:
            new_elem_num = newMesh->AddHex(indices, n + 1); break;
         case ET_PENTA6: 
         // case ET_PENTA15: // 2nd order wedge TODO
            new_elem_num = newMesh->AddWedge(indices, n + 1); break;
         case ET_PYRA5: 
         // case ET_PYRA13: // 2nd order pyramid TODO
            new_elem_num = newMesh->AddPyramid(indices, n + 1);
            break;
			default:
				assert(false);
			}

         // newMesh->GetElementCenter(new_elem_num, elem_center);
         // elem_center_std[0] = elem_center(0);
         // elem_center_std[1] = elem_center(1);
         // elem_center_std[2] = elem_center(2);
         // newMesh->SetFEBioElementID(elem_center_std, el.GetID());
		}
	}

	for (int n = 0; n < febioMesh.FacetSets(); ++n)
	{
		FEFacetSet& surf = febioMesh.FacetSet(n);
		int faces = surf.Faces();
		newMesh->AssignSurfaceAttribute(&surf, n + 1);
		for (int i = 0; i < faces; ++i)
		{
         std::fill(indices, indices + FEElement::MAX_NODES, 0); // reset indices to all 0
			FEFacetSet::FACET& face = surf.Face(i);
			for (int j = 0; j < face.Nodes(); j++) {
				indices[j] = face.node[j];
			}

         switch (face.ntype)
         {
         case FEFacetSet::FACET::TRI3: 
         case FEFacetSet::FACET::TRI6:
            newMesh->AddBdrTriangle(indices, n + 1); break;
         case FEFacetSet::FACET::QUAD4:
         case FEFacetSet::FACET::QUAD9:
            newMesh->AddBdrQuad(indices, n + 1); break;
         default:
            assert(false);
         }
		}
	}

   newMesh->FinalizeTopology(false);
	newMesh->Finalize(false, true);

	return newMesh;
}