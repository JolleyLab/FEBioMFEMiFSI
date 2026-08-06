#include "MFEMParImport.h"

// FEBio Headers
#include <FECore/XMLReader.h>
#include <FEBioXML/FEBioIncludeSection.h>
#include <FEBioXML/FEBioModuleSection.h>
#include <FEBioXML/FEBioControlSection4.h>
#include <FEBioXML/FEBioGlobalsSection.h>
#include <FEBioXML/FEBioMaterialSection.h>
#include <FEBioXML/FEBioGeometrySection.h>
#include <FEBioXML/FEBioBoundarySection3.h>
#include <FEBioXML/FEBioLoadsSection.h>
#include <FEBioXML/FEBioContactSection.h>
#include <FEBioXML/FEBioConstraintsSection.h>
#include <FEBioXML/FEBioInitialSection3.h>
#include <FEBioXML/FEBioLoadDataSection.h>
#include <FEBioXML/FEBioOutputSection.h>
#include <FEBioXML/FEBioStepSection4.h>
#include <FEBioXML/FEBioDiscreteSection.h>
#include <FEBioXML/FEBioMeshDataSection.h>
#include <FEBioXML/FEBioCodeSection.h>
#include <FEBioXML/FEBioRigidSection4.h>
#include <FEBioXML/FEBioMeshAdaptorSection.h>
#include <FEBioXML/FEBioMeshSection4.h>
#include <FEBioXML/FEBioMeshDomainsSection4.h>
#include <FEBioXML/FEBioStepSection3.h>

// MFEM Headers
#include "MFEMBoundarySection.h"
#include "MFEMLoadSection.h"
#include "MFEMMeshSection.h"
#include "MFEMMeshDomainsSection.h"

MFEMParImport::MFEMParImport() : FEBioImport()
{
	// Constructor implementation
}

MFEMParImport::~MFEMParImport()
{
   // Destructor implementation
   FEBioImport::~FEBioImport();
}

bool MFEMParImport::Load(FEModel& fem, const char* szfile)
{
   if (m_builder == nullptr)
	{
		m_builder = new FEModelBuilder(fem);
	}

	// intialize some variables
	m_szdmp[0] = 0;
	m_szlog[0] = 0;
	m_szplt[0] = 0;

	m_data.clear();

	// extract the path
	strcpy(m_szpath, szfile);
	char* ch = strrchr(m_szpath, '\\');
	if (ch==0) ch = strrchr(m_szpath, '/');
	if (ch==0) m_szpath[0] = 0; else *(ch+1)=0;

	// clean up
	fem.GetMesh().ClearDataMaps();

	// read the file
	if (ReadFile(szfile) == false) return false;

	// finish building
	try {
		bool b = m_builder->Finish();
		if (b == false) return errf("FAILED building FEBio model.");
	}
	catch (std::exception e)
	{
		const char* szerr = e.what();
		if (szerr == nullptr) szerr = "(unknown exception)";
		return errf("%s", e.what());
	}	
	catch (...)
	{
		return errf("unknown exception.");
	}

	return true;
}

bool MFEMParImport::ReadFile(const char* szfile, bool broot)
{
   // Open the XML file
	XMLReader xml;
	if (xml.Open(szfile) == false) return errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);

	// Find the root element
	XMLTag tag;
	try
	{
		if (xml.FindTag("febio_spec", tag) == false) return errf("FATAL ERROR: febio_spec tag was not found. This is not a valid input file.\n\n");
	}
	catch (...)
	{
		return errf("An error occured while finding the febio_spec tag.\nIs this a valid FEBio input file?\n\n");
	}

	// parse the file
	try
	{
		// get the version number
		ParseVersion(tag);

		// FEBio4 only supports file version 1.2, 2.0, 2.5, 3.0, and 4.0
		int nversion = GetFileVersion();
      if (nversion != 0x0400) throw InvalidVersion(); // assumes 4.0

		// build the file section map based on the version number
		// BuildFileSectionMap(nversion);
      BuildMFEMFileSectionMap();

		// parse the file
		++tag;

		do
		{
			// try to find a section parser
			FEFileSectionMap::iterator is = m_map.find(tag.Name());

			// make sure we found a section reader
			if (is == m_map.end()) throw XMLReader::InvalidTag(tag);

			// parse the section
         is->second->Parse(tag);

			// go to the next tag
			++tag;
		}
		while (!tag.isend());
	}
	// --- XML Reader Exceptions ---
	catch (XMLReader::Error& e)
	{
		return errf("%s", e.what());
	}
	// --- FEBioImport Exceptions ---
	catch (FEFileException& e)
	{
		return errf("%s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
	}
	// --- Exception from DataStore ---
	catch (UnknownDataField& e)
	{
		return errf("\"%s\" is not a valid field variable name (line %d)\n", e.what(), xml.GetCurrentLine()-1);
	}
	// std::exception
	catch (std::exception e)
	{
		const char* szerr = e.what();
		if (szerr == nullptr) szerr = "(unknown exception)";
		return errf("%s", szerr);
	}
	// --- Unknown exceptions ---
	catch (...)
	{
		return errf("unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	// close the XML file
	xml.Close();

	// we're done!
	return true;
}

void MFEMParImport::BuildMFEMFileSectionMap()
{
   // do not allow unknown attributes
   SetStopOnUnknownAttribute(true);

   m_map["Control"    ] = new FEBioControlSection4     (this);
   m_map["Material"   ] = new FEBioMaterialSection3    (this);

   // m_map["Mesh"       ] = new FEBioMeshSection4        (this);
   m_map["Mesh"       ] = new MFEMMeshSection          (this);

   // m_map["MeshDomains"] = new FEBioMeshDomainsSection4 (this);
   m_map["MeshDomains"] = new MFEMMeshDomainsSection   (this, solid_mesh);

   m_map["Include"    ] = new FEBioIncludeSection      (this);
   m_map["Initial"    ] = new FEBioInitialSection3     (this);

   // m_map["Boundary"   ] = new FEBioBoundarySection3    (this);
   m_map["Boundary"   ] = new MFEMBoundarySection      (this);

   // m_map["Loads"      ] = new FEBioLoadsSection3       (this);
   m_map["Loads"      ] = new MFEMLoadSection          (this);

   m_map["Contact"    ] = new FEBioContactSection4     (this);
   m_map["Discrete"   ] = new FEBioDiscreteSection25   (this);
   m_map["Constraints"] = new FEBioConstraintsSection25(this);
   m_map["Code"       ] = new FEBioCodeSection         (this); // added in FEBio 2.4 (experimental feature!)
   m_map["MeshData"   ] = new FEBioMeshDataSection4    (this);	// added in febio4
   m_map["LoadData"   ] = new FEBioLoadDataSection3    (this);
   m_map["Rigid"      ] = new FEBioRigidSection4       (this); // added in FEBio 4.0
   m_map["Step"       ] = new FEBioStepSection4        (this);
   m_map["MeshAdaptor"] = new FEBioMeshAdaptorSection  (this);	// added in FEBio 3.0
}
