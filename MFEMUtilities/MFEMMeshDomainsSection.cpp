#include "MFEMMeshDomainsSection.h"

#include <FECore/FEDomain.h>
#include <FECore/FEMaterial.h>


void MFEMMeshDomainsSection::Parse(XMLTag& tag)
{
	// build the node ID lookup table
	BuildNLT();
	
	// read all sections
	if (tag.isleaf() == false)
	{
		++tag;
		do
		{
			if      (tag == "SolidDomain") ParseSolidDomainSection(tag);
			else if (tag == "ShellDomain") ParseShellDomainSection(tag);
			else if (tag == "BeamDomain" ) ParseBeamDomainSection(tag);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
	}

	// let's build the part
	FEBModel& feb = GetBuilder()->GetFEBModel();
	FEBModel::Part* part = feb.GetPart(0); assert(part);
	if (feb.BuildPart(*GetFEModel(), *part, false) == false)
	{
		throw XMLReader::Error("Failed building parts.");
	}

	// tell the file reader to rebuild the node ID table
	GetBuilder()->BuildNodeList();

	// At this point the mesh is completely read in.
	// allocate material point data
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		dom.CreateMaterialPointData();
	}

	// Now we can allocate the degrees of freedom.
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	fem.GetMesh().SetDOFS(MAX_DOFS);
}

void MFEMMeshDomainsSection::BuildNLT()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	FEBModel& feb = GetBuilder()->GetFEBModel();
	FEBModel::Part* part = feb.GetPart(0); assert(part);
	if (part == nullptr) return;

	// build node-index lookup table
	int noff = -1, maxID = 0;
	int N0 = mesh.Nodes();
	int NN = part->Nodes();
	for (int i = 0; i < NN; ++i)
	{
		int nid = part->GetNode(i).id;
		if ((noff < 0) || (nid < noff)) noff = nid;
		if (nid > maxID) maxID = nid;
	}
	m_NLT.assign(maxID - noff + 1, -1);
	for (int i = 0; i < NN; ++i)
	{
		int nid = part->GetNode(i).id - noff;
		m_NLT[nid] = i + N0;
	}

	m_noff = noff;
}

void MFEMMeshDomainsSection::ParseSolidDomainSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FEBModel& feb = GetBuilder()->GetFEBModel();
	FEBModel::Part* part = feb.GetPart(0); assert(part);
	if (part == nullptr) throw XMLReader::InvalidTag(tag);

	// get the domain name
	const char* szname = tag.AttributeValue("name");
	FEBModel::Domain* partDomain = part->FindDomain(szname);
	if (partDomain == nullptr) throw XMLReader::InvalidAttributeValue(tag, "name", szname);

	// get the material name
	const char* szmat = tag.AttributeValue("mat");
	FEMaterial* mat = fem.FindMaterial(szmat);
	if (mat == nullptr) throw XMLReader::InvalidAttributeValue(tag, "mat", szmat);

	// set the material name
	partDomain->SetMaterialName(szmat);

	// get the (optional) type attribute
	const char* sztype = tag.AttributeValue("type", true);

	// --- build the domain --- 
	// we'll need the kernel for creating domains
	FECoreKernel& febio = FECoreKernel::GetInstance();

	// element count
	int elems = partDomain->Elements();

	// get the element spect
	FE_Element_Spec spec = partDomain->ElementSpec();
   spec.eclass = FE_ELEM_SOLID;

	// create the domain
	FEDomain* dom = nullptr;
	if (sztype)
	{
		// if the type attribute is defined, try to allocate the domain class directly. 
		dom = febio.CreateDomainExplicit(FESOLIDDOMAIN_ID, sztype, &fem);
		if (dom == nullptr) throw XMLReader::InvalidAttributeValue(tag, sztype);
		dom->SetMaterial(mat);
	}
	else
	{
		// if not, then use "old" logic, which tries to match the domain using the 
		// domain factories. 
		dom = febio.CreateDomain(spec, &mesh, mat);
		if (dom == 0) throw XMLReader::InvalidTag(tag);
	}

	// add it to the mesh
	mesh.AddDomain(dom);

	// Allocate elements
	if (dom->Create(elems, spec) == false)
	{
		throw XMLReader::InvalidTag(tag);
	}

	// assign the material
	dom->SetMatID(mat->GetID() - 1);

	// get the part name
	string partName = part->Name();
	if (partName.empty() == false) partName += ".";

	string domName = partName + partDomain->Name();
	dom->SetName(domName);

	// read additional parameters
	if (tag.isleaf() == false)
	{
		ReadParameterList(tag, dom);
	}

	// process element data
	for (int j = 0; j < elems; ++j)
	{
		const FEBModel::ELEMENT& domElement = partDomain->GetElement(j);

		FEElement& el = dom->ElementRef(j);
		el.SetID(domElement.id);

      // set element type
      int mfem_elem_type = solid_mesh->GetElementType(domElement.id-1); // MFEM_ID = FEBIO_ID - 1
      // int nnodes = domElement.node
      SetElementSpec(mfem_elem_type,spec);
      el.SetType(spec.etype);

		int ne = el.Nodes();
		for (int n = 0; n < ne; ++n) el.m_node[n] = m_NLT[domElement.node[n] - m_noff];
	}
}

void MFEMMeshDomainsSection::SetElementSpec(int mfem_type, FE_Element_Spec &spec)
{
   // Assumes linear elements
   if (mfem_type == mfem::Element::TETRAHEDRON)
   {
      spec.eshape = FE_Element_Shape::ET_TET4;
      spec.etype = FE_Element_Type::FE_TET4G4;
   }
   else if (mfem_type == mfem::Element::HEXAHEDRON)
   {
      spec.eshape = FE_Element_Shape::ET_HEX8;
      spec.etype = FE_Element_Type::FE_HEX8G8;
   }
   else if (mfem_type == mfem::Element::WEDGE)
   {
      spec.eshape = FE_Element_Shape::ET_PENTA6;
      spec.etype = FE_Element_Type::FE_PENTA6G6;
   }
   else if (mfem_type == mfem::Element::PYRAMID)
   {
      spec.eshape = FE_Element_Shape::ET_PYRA5;
      spec.etype = FE_Element_Type::FE_PYRA5G8;
   }
   else
   {
      MFEM_ABORT("Unsupported element type: " << mfem_type);
   }
}

void MFEMMeshDomainsSection::ParseShellDomainSection(XMLTag& tag)
{
   MFEM_ABORT("TODO");
}

void MFEMMeshDomainsSection::ParseBeamDomainSection(XMLTag& tag)
{
   MFEM_ABORT("TODO");
}