#pragma once

#include <FEBioXML/FileImport.h>

#include "mfem.hpp"

class MFEMMeshDomainsSection : public FEFileSection
{
public:
   MFEMMeshDomainsSection(FEFileImport* imp, mfem::Mesh* solid_mesh) : FEFileSection(imp), solid_mesh(solid_mesh) { m_noff = 0; };

   void Parse(XMLTag& tag);

protected:
	void ParseSolidDomainSection(XMLTag& tag);
   void ParseShellDomainSection(XMLTag& tag);
	void ParseBeamDomainSection(XMLTag& tag);
private:
	void BuildNLT();
   void SetElementSpec(int mfem_type, FE_Element_Spec &spec);

private:
	std::vector<int>	m_NLT;
	int					m_noff;
	mfem::Mesh*			solid_mesh;
};