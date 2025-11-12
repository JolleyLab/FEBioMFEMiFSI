#pragma once

#include <FEBioXML/FileImport.h>

class MFEMMeshSection : public FEFileSection
{
public:
   MFEMMeshSection(FEFileImport* imp) : FEFileSection(imp) { m_maxNodeId = 0; };

   void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag, FEBModel::Part* part);
	void ParseSurfaceSection    (XMLTag& tag, FEBModel::Part* part);
	void ParseElementSection    (XMLTag& tag, FEBModel::Part* part);
private:
   int m_maxNodeId;
};