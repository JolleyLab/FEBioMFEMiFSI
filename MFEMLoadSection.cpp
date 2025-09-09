#include "MFEMLoadSection.h"

#include <FECore/FENodalBC.h>

#include "mfem.hpp"
#include "MFEMLoad.h"


void MFEMLoadSection::Parse(XMLTag& tag)
{
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "nodal_load"  ) ParseMFEMNodalLoad  (tag);
		else if (tag == "surface_load") ParseMFEMSurfaceLoad(tag);
		else if (tag == "body_load"   ) ParseMFEMBodyLoad   (tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

void MFEMLoadSection::ParseMFEMBodyLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
	const char* sztype = tag.AttributeValue("type");
	FEBodyLoad* pbl = fecore_new<FEBodyLoad>(sztype, &fem);
	if (pbl == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// see if a name was defined
	const char* szname = tag.AttributeValue("name", true);
	if (szname) pbl->SetName(szname);

	// see if a specific domain was referenced
	const char* szpart = tag.AttributeValue("elem_set", true);
	// if (szpart)
	// {
	// 	FEMesh& mesh = fem.GetMesh();
	// 	FEElementSet* elset = mesh.FindElementSet(szpart);
	// 	if (elset == 0) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szpart);
	// 	pbl->SetDomainList(elset);
	// }

	// read the mfem domain attribute
	MFEMBodyLoad* mfem_bl = dynamic_cast<MFEMBodyLoad*>(pbl);
	MFEM_ASSERT(mfem_bl,"Error: not a MFEM body load");
	int mfem_dom_attr;
	tag.AttributeValue("mfem_dom_attr",mfem_dom_attr,true);
	mfem_bl->SetDomainAttr(mfem_dom_attr);

	ReadParameterList(tag, pbl);
	GetBuilder()->AddModelLoad(pbl);
}

void MFEMLoadSection::ParseMFEMNodalLoad(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the type
	const char* sztype = tag.AttributeValue("type");

	// create nodal load
	FENodalLoad* pfc = fecore_new<FENodalLoad>(sztype, GetFEModel());
	if (pfc == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type");

	// read (optional) name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) pfc->SetName(szname);

	// get the node set
	const char* szset = tag.AttributeValue("node_set");
	// FENodeSet* nodeSet = GetBuilder()->FindNodeSet(szset);
	// if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// // assign the node set
	// pfc->SetNodeSet(nodeSet);

	// add it to the model
	GetBuilder()->AddNodalLoad(pfc);

	// read parameters
	ReadParameterList(tag, pfc);
}

void MFEMLoadSection::ParseMFEMSurfaceLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create surface load
	string typeString = tag.AttributeValue("type");
	FESurfaceLoad* psl = fecore_new<FESurfaceLoad>(typeString.c_str(), &fem);
	if (psl == 0) throw XMLReader::InvalidTag(tag);

	// read (optional) name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) psl->SetName(szname);

	// read required surface attribute
	const char* surfaceName = tag.AttributeValue("surface");
	// FEFacetSet* pface = mesh.FindFacetSet(surfaceName);
	// if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", surfaceName);
   
   // create dummy facet set
   FEFacetSet* pface = new FEFacetSet(0); // with 0 facets
   mesh.AddFacetSet(pface); // pface is deleted by the mesh

	// create a surface from this facet set
	FESurface* psurf = fecore_alloc(FESurface, &fem);
	GetBuilder()->BuildSurface(*psurf, *pface);

	// assign it
	// mesh.AddSurface(psurf);
	psl->SetSurface(psurf);

   // read the mfem bdr attribute
   MFEMSurfaceLoad* mfem_sl = dynamic_cast<MFEMSurfaceLoad*>(psl);
   MFEM_ASSERT(mfem_sl,"Error: not a MFEM surface load");
   int mfem_bdr_attr;
   tag.AttributeValue("mfem_bdr_attr",mfem_bdr_attr,true);
   mfem_sl->SetBoundaryID(mfem_bdr_attr);

	// add it to the model
	GetBuilder()->AddSurfaceLoad(psl);

	// read the parameters
	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if (ReadParameter(tag, psl) == false)
			{
				// some special handling for fluidnormalvelocity
				if ((tag == "value") && (strcmp(typeString.c_str(), "fluid normal velocity") == 0))
				{
					// The value parameter no longer exists, but for backward compatibility
					// we map it to the "velocity" parameter.

					// get the surface data attribute
					const char* szsurfdata = tag.AttributeValue("surface_data");

					// find the surface map
					FEDataMap* surfMap = mesh.FindDataMap(szsurfdata);
					if (surfMap == nullptr) throw XMLReader::InvalidAttributeValue(tag, "surface_data");

					// get the velocity parameter
					FEParam* p = psl->GetParameter("velocity");
					FEParamDouble& v = p->value<FEParamDouble>();

					// we'll map the velocity value to the map's scale factor
					double s = 1.0;
					if (v.isConst()) s = v.constValue(); else throw XMLReader::InvalidTag(tag);

					// create map and assign
					FEMappedValue* map = new FEMappedValue(&fem);
					map->setDataMap(surfMap);
					map->setScaleFactor(s);
					v.setValuator(map);
				}
				else throw XMLReader::InvalidTag(tag);
			}
			++tag;
		} 
		while (!tag.isend());
	}
}