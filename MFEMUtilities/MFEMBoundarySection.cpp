#include "MFEMBoundarySection.h"

#include <FECore/FENodalBC.h>

#include "mfem.hpp"
#include "../MFEMCore/MFEMBoundaryCondition.h"


void MFEMBoundarySection::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if (tag == "bc") ParseMFEMBC(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

void MFEMBoundarySection::ParseMFEMBC(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();

	// get the type string
	const char* sztype = tag.AttributeValue("type");

	// handle "rigid" bc separately
	if (strcmp(sztype, "rigid") == 0)
	{
      MFEM_ABORT("TODO");
		// ParseBCRigid(tag);
		return;
	}

	// handle linear constraints
	if (strcmp(sztype, "linear constraint") == 0)
	{
      MFEM_ABORT("TODO");
		// ParseLinearConstraint(tag);
		return;
	}

	// create the boundary condition
	FEBoundaryCondition* pbc = fecore_new<FEBoundaryCondition>(sztype, fem);
	if (pbc == 0) throw XMLReader::InvalidAttributeValue(tag, "type");

	// read the (optional) name 
	const char* szname = tag.AttributeValue("name", true);
	if (szname) pbc->SetName(szname);

   // read the mfem bdr attribute
   MFEMBoundaryCondition* mfem_bc = dynamic_cast<MFEMBoundaryCondition*>(pbc);
   MFEM_ASSERT(mfem_bc,"Error: not a MFEM boundary condition");
   int mfem_bdr_attr;
   tag.AttributeValue("mfem_bdr_attr",mfem_bdr_attr,true);
   mfem_bc->SetBoundaryID(mfem_bdr_attr);

	// get the node set
	if (dynamic_cast<FENodalBC*>(pbc))
	{
		FENodalBC* nbc = dynamic_cast<FENodalBC*>(pbc); assert(nbc);

		// read required node_set attribute
		const char* szset = tag.AttributeValue("node_set");
		// FENodeSet* nodeSet = GetBuilder()->FindNodeSet(szset);
		// if (nodeSet == nullptr) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

		// nbc->SetNodeSet(nodeSet);

      MFEM_ABORT("TODO");
	}

	// get the surface
	if (dynamic_cast<FESurfaceBC*>(pbc))
	{
		FESurfaceBC* sbc = dynamic_cast<FESurfaceBC*>(pbc);

		// read required surface attribute
		const char* surfaceName = tag.AttributeValue("surface");
		// FEFacetSet* pface = mesh.FindFacetSet(surfaceName);
		// if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", surfaceName);

      // // create dummy facet set
      // FEFacetSet* pface = new FEFacetSet(0); // with 0 facets
      // mesh.AddFacetSet(pface); // pface is deleted by the mesh

		// // create a surface from this facet set
		// FESurface* psurf = fecore_alloc(FESurface, fem);
		// GetBuilder()->BuildSurface(*psurf, *pface);

		// // assign it
		// mesh.AddSurface(psurf);
		// sbc->SetSurface(psurf);

      mfem_bc->SetSurfaceName(surfaceName);
	}

	// add this boundary condition to the current step
	GetBuilder()->AddBC(pbc);

	// Read the parameter list
	ReadParameterList(tag, pbc);
}