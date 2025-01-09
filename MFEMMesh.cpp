#include "MFEMMesh.h"
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include "mfem.hpp"
using namespace mfem;

// This function constructs a mfem mesh from the FEBioMesh.
// Note that this currently assumes a 2D geometry! 
// TODO: Generalize this for 3D meshes.
mfem::Mesh* CreateMFEMMesh(FEMesh& febioMesh)
{
	int nodes = febioMesh.Nodes();
	int elems = febioMesh.Elements();

	// In FEBio, it is not necessary to provide the entire boundary. 
	// I don't know if this will be an issue for mfem. 
	int bdrElems = 0;
	for (int i = 0; i < febioMesh.SegmentSets(); ++i)
	{
		bdrElems += febioMesh.SegmentSet(i).Segments();
	}

	Mesh* mesh = new Mesh(2, nodes, elems, bdrElems, 2);
	for (int i = 0; i < nodes; ++i)
	{
		FENode& node = febioMesh.Node(i);
		vec3d& pos = node.m_r0;
		double coords[3] = { pos.x, pos.y, pos.z };
		mesh->AddVertex(coords);
	}

	int indices[FEElement::MAX_NODES] = { 0 };
	for (int i = 0; i < febioMesh.SegmentSets(); ++i)
	{
		FESegmentSet& set = febioMesh.SegmentSet(i);
		for (int j = 0; j < set.Segments(); ++j)
		{
			FESegmentSet::SEGMENT& line = set.Segment(j);
			indices[0] = line.node[0];
			indices[1] = line.node[1];
			mesh->AddBdrSegment(indices[0], indices[1], i + 1);
		}
	}

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
			case ET_QUAD4: mesh->AddQuad(indices[0], indices[1], indices[2], indices[3], n + 1); break;
			default:
				assert(false);
			}
		}
	}

	mesh->FinalizeTopology(false);
	mesh->Finalize(false, true);

	return mesh;
}
