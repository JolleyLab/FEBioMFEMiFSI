#pragma once
#include <FECore/FECoreFactory.h>
#include <FECore/FEShellDomain.h>
#include "MFEMMaterial.h"

class MFEMDomain : public FEShellDomainOld
{
public:
	MFEMDomain(FEModel* fem) : FEShellDomainOld(fem), dofs(fem) {}

	const FEDofList& GetDOFList() const override
	{
		return dofs;
	}

	void SetMaterial(FEMaterial* pmat) override
	{
		mat = dynamic_cast<MFEMMaterial*>(pmat);
		FEShellDomainOld::SetMaterial(pmat);
	}

	FEMaterial* GetMaterial() override { return mat; }

private:
	FEDofList dofs;
	MFEMMaterial* mat;
};

class MFEMDomainFactory : public FEDomainFactory
{
public:
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
