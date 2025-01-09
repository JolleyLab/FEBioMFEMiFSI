#include "MFEMDomain.h"
#include <FECore/FEModel.h>
#include "MFEMMaterial.h"

FEDomain* MFEMDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel& fem = *pmat->GetFEModel();
	FE_Element_Class eclass = spec.eclass;
	const char* sztype = 0;
	if (dynamic_cast<MFEMMaterial*>(pmat))
	{
		if (eclass == FE_ELEM_SHELL) sztype = "mfem domain";
	}

	if (sztype)
	{
		FEDomain* pd = fecore_new<FEShellDomainOld>(sztype, &fem);
		if (pd) pd->SetMaterial(pmat);
		return pd;
	}
	else return nullptr;
}
