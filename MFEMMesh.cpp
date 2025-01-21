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
	FESurface* boundaryMesh = febioMesh.ElementBoundarySurface();

	int nodes = febioMesh.Nodes();
	int elems = febioMesh.Elements();
	int nbndrElems = boundaryMesh->Elements();
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
		for (int i = 0; i < elemsInDomain; i++) {
			FEElement& el = dom.ElementRef(i);
			for (int j = 0; j < el.Nodes(); j++) {
				indices[j] = el.m_node[j];
			}

			switch (el.Shape())
			{
			case ET_HEX8: newMesh->AddHex(indices, n + 1); break;
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

			if (face.ntype == FEFacetSet::FACET::QUAD4)
			{
				newMesh->AddBdrQuad(indices, n + 1);
			}
		}
	}

	// newMesh->FinalizeTopology(false);
   newMesh->FinalizeTopology();
	newMesh->Finalize(false, true);

	return newMesh;
}
