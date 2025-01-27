#include "MFEMMesh.h"
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESurface.h>
#include "mfem.hpp"
#include <math.h>
using namespace mfem;

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
	for (int n = 0; n < febioMesh.Domains(); ++n)
	{
		FEDomain& dom = febioMesh.Domain(n);
		int elemsInDomain = dom.Elements();
      newMesh->AssignDomainAttribute(&dom, n + 1);
		for (int i = 0; i < elemsInDomain; i++) {
			FEElement& el = dom.ElementRef(i);
			for (int j = 0; j < el.Nodes(); j++) {
				indices[j] = el.m_node[j];
			}

			switch (el.Shape())
			{
         case ET_TET4: 
         case ET_TET10:
            newMesh->AddTet(indices, n + 1); break;
         case ET_HEX8: 
         case ET_HEX27:
            newMesh->AddHex(indices, n + 1); break;
         case ET_PENTA6: 
         // case ET_PENTA15: // 2nd order wedge TODO
            newMesh->AddWedge(indices, n + 1); break;
         case ET_PYRA5: 
         // case ET_PYRA13: // 2nd order pyramid TODO
            newMesh->AddPyramid(indices, n + 1);
            assert(false); // no routines for computing d2Ndxx in MFEM yet
            break;
			default:
				assert(false);
			}
		}
	}

	for (int n = 0; n < febioMesh.FacetSets(); ++n)
	{
		FEFacetSet& surf = febioMesh.FacetSet(n);
		int faces = surf.Faces();
		newMesh->AssignSurfaceAttribute(&surf, n + 1);
		for (int i = 0; i < faces; ++i)
		{
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
