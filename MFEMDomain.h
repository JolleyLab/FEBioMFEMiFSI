#pragma once
#include <FECore/FECoreFactory.h>
#include <FECore/FEShellDomain.h>
#include <FECore/FESolidDomain.h>
#include "MFEMMaterial.h"

class MFEMDomain2D : public FEShellDomainOld
{
public:
	MFEMDomain2D(FEModel* fem) : FEShellDomainOld(fem), dofs(fem), mat(nullptr) {}

	const FEDofList& GetDOFList() const override
	{
		return dofs;
	}

	void SetMaterial(FEMaterial* pmat) override
	{
      if (dynamic_cast<MFEMMaterial*>(pmat)) 
      {
         mat = dynamic_cast<MFEMMaterial*>(pmat);
         FEShellDomainOld::SetMaterial(pmat);
      }
      else if (dynamic_cast<MFEMUncoupledFEBioMaterial*>(pmat)) 
      {
         mat = dynamic_cast<MFEMUncoupledFEBioMaterial*>(pmat);
         FEShellDomainOld::SetMaterial(pmat);
      }
      else assert(false);
	}

	FEMaterial* GetMaterial() override { return mat; }

private:
	FEDofList dofs;
	FEMaterial* mat;
};


class MFEMSolidDomain : public FESolidDomain
{
public:
	MFEMSolidDomain(FEModel* fem) : FESolidDomain(fem), dofs(fem), mat(nullptr) {}

	const FEDofList& GetDOFList() const override
	{
		return dofs;
	}

	void SetMaterial(FEMaterial* pmat) override
	{
      if (dynamic_cast<MFEMMaterial*>(pmat)) 
      {
         mat = dynamic_cast<MFEMMaterial*>(pmat);
         FESolidDomain::SetMaterial(pmat);
      }
      else if (dynamic_cast<MFEMUncoupledFEBioMaterial*>(pmat)) 
      {
         mat = dynamic_cast<MFEMUncoupledFEBioMaterial*>(pmat);
         FESolidDomain::SetMaterial(pmat);
      }
      else assert(false);
	}

	FEMaterial* GetMaterial() override { return mat; }

private:
	FEDofList dofs;
	FEMaterial* mat;
};

class MFEMDomainFactory : public FEDomainFactory
{
public:
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
