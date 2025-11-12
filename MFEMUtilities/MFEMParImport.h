#pragma once

#include <FEBioXML/FEBioImport.h>

#include "mfem.hpp"

class MFEMParImport : public FEBioImport
{
public:
	// constructor
   MFEMParImport();

   // destructor
   virtual ~MFEMParImport();

   // open the file
   virtual bool Load(FEModel& fem, const char* szfile);

   // read file contents
   virtual bool ReadFile(const char* szfile, bool broot = false);

   void SetSolidMesh(mfem::Mesh* mesh) { solid_mesh = mesh; }

private:
   void BuildMFEMFileSectionMap();

   mfem::Mesh *solid_mesh = nullptr;
};