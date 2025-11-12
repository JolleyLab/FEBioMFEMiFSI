#pragma once

#include <FEBioXML/FileImport.h>

class MFEMLoadSection : public FEFileSection
{
public:
   MFEMLoadSection(FEFileImport* imp) : FEFileSection(imp) {};

   void Parse(XMLTag& tag);
	void ParseMFEMNodalLoad(XMLTag& tag);
   void ParseMFEMSurfaceLoad(XMLTag& tag);
   void ParseMFEMBodyLoad(XMLTag& tag);
};