#pragma once

#include <FEBioXML/FEBioBoundarySection.h>

class MFEMBoundarySection : public FEBioBoundarySection
{
public:
   MFEMBoundarySection(FEFileImport* imp) : FEBioBoundarySection(imp) {};

   void Parse(XMLTag& tag);
	void ParseMFEMBC(XMLTag& tag);
};