#include "MFEMMeshSection.h"

#include "mfem.hpp"


void MFEMMeshSection::Parse(XMLTag& tag)
{
	FEModelBuilder* builder = GetBuilder();
	builder->m_maxid = 0;

	// create a default part
	// NOTE: Do not specify a name for the part, otherwise
	//       all lists will be given the name: partname.listname
	FEBModel& feb = builder->GetFEBModel();
	assert(feb.Parts() == 0);
	FEBModel::Part* part = feb.AddPart("");

	// read all sections
	++tag;
	do
	{
		if      (tag == "Nodes"      ) ParseNodeSection       (tag, part);
		else if (tag == "Elements"   ) ParseElementSection    (tag, part);
		else if (tag == "Surface"    ) ParseSurfaceSection    (tag, part);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

void MFEMMeshSection::ParseNodeSection(XMLTag& tag, FEBModel::Part* part)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// see if this list defines a set
	const char* szname = tag.AttributeValue("name", true);
	FEBModel::NodeSet* ps = 0;
	if (szname)
	{
		ps = new FEBModel::NodeSet(szname);
		part->AddNodeSet(ps);
	}

	// allocate node

	vector<FEBModel::NODE> node; node.reserve(10000);
	vector<int> nodeList; nodeList.reserve(10000);

	// read nodal coordinates
	++tag;
	do {
		// nodal coordinates
		FEBModel::NODE nd;
		value(tag, nd.r);

		// get the nodal ID
		tag.AttributeValue("id", nd.id);

		// make sure node IDs are incrementing
		// if (nd.id <= m_maxNodeId) throw XMLReader::InvalidAttributeValue(tag, "id");
		// m_maxNodeId = nd.id;

		// add it to the pile
		node.push_back(nd);
		nodeList.push_back(nd.id);

		// go on to the next node
		++tag;
	} while (!tag.isend());

	// add nodes to the part
	part->AddNodes(node);

	// If a node set is defined add these nodes to the node-set
	if (ps) ps->SetNodeList(nodeList);
}

void MFEMMeshSection::ParseElementSection(XMLTag& tag, FEBModel::Part* part)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the (required!) name
	const char* szname = tag.AttributeValue("name");

	// // get the element spec
	// const char* sztype = tag.AttributeValue("type");
	// FE_Element_Spec espec = GetBuilder()->ElementSpec(sztype);
	// if (FEElementLibrary::IsValid(espec) == false) throw FEBioImport::InvalidElementType();

	// make sure the domain does not exist yet
	FEBModel::Domain* dom = part->FindDomain(szname);
	if (dom)
	{
		stringstream ss;
		ss << "Duplicate part name found : " << szname;
		throw std::runtime_error(ss.str());
	}

	// create the new domain
	// dom = new FEBModel::Domain(espec);
   dom = new FEBModel::Domain();
	if (szname) dom->SetName(szname);

	// add domain it to the mesh
	part->AddDomain(dom);

	// for named domains, we'll also create an element set
	FEBModel::ElementSet* pg = 0;
	if (szname)
	{
		pg = new FEBModel::ElementSet(szname);
		part->AddElementSet(pg);
	}

	dom->Reserve(10000);
	vector<int> elemList; elemList.reserve(10000);

	// keep track of largest ID
	// we need to enforce that element IDs are increasing
	// (This is currently only done for each domain. Need to modify this
	// so it's done on the whole model.)
	int maxID = -1;

	// read element data
	++tag;
	do
	{
		FEBModel::ELEMENT el;

		// get the element ID
		tag.AttributeValue("id", el.id);

		if ((maxID == -1) || (el.id > maxID)) maxID = el.id;
		else throw XMLReader::InvalidAttributeValue(tag, "id");

		// read the element data
		tag.value(el.node, FEElement::MAX_NODES);

		dom->AddElement(el);
		elemList.push_back(el.id);

		// go to next tag
		++tag;
	} while (!tag.isend());

	// set the element list
	if (pg) pg->SetElementList(elemList);
}

void MFEMMeshSection::ParseSurfaceSection(XMLTag& tag, FEBModel::Part* part)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// // see if this surface was already defined
	FEBModel::Surface* ps = part->FindSurface(szname);
	// if (ps) throw FEBioImport::RepeatedSurface(szname);

	// allocate storage for faces
	ps = new FEBModel::Surface(szname);
	part->AddSurface(ps);

	if (tag.isleaf() || tag.isempty()) return;

	// count nr of faces
	int faces = tag.children();
	ps->Create(faces);

	// read faces
	++tag;
	for (int i = 0; i < faces; ++i)
	{
		FEBModel::FACET& face = ps->GetFacet(i);

		// get the ID (although we don't really use this)
		tag.AttributeValue("id", face.id);

		// set the facet type
		if (tag == "quad4") face.ntype = 4;
		else if (tag == "tri3") face.ntype = 3;
		else if (tag == "tri6") face.ntype = 6;
		else if (tag == "tri7") face.ntype = 7;
		else if (tag == "quad8") face.ntype = 8;
		else if (tag == "quad9") face.ntype = 9;
		else throw XMLReader::InvalidTag(tag);

		// we assume that the facet type also defines the number of nodes
		int N = face.ntype;
		tag.value(face.node, N);

		++tag;
	}
}