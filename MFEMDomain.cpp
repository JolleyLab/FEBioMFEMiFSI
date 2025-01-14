#include "MFEMDomain.h"
#include <FECore/FEModel.h>
#include "MFEMMaterial.h"

FEDomain* MFEMDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel& fem = *pmat->GetFEModel();
	FE_Element_Class eclass = spec.eclass;
	FEDomain* dom = nullptr;
	if (dynamic_cast<MFEMMaterial*>(pmat))
	{
		if (eclass == FE_ELEM_SHELL)
		{
			dom = fecore_new<FEShellDomainOld>("mfem domain 2d", &fem);
		}
		else if (eclass == FE_ELEM_SOLID)
		{
			dom = fecore_new<FESolidDomain>("mfem domain 3d", &fem);
		}
		if (dom) dom->SetMaterial(pmat);
	}
	return dom;
}
