#include "MFEMDomain.h"
#include <FECore/FEModel.h>
#include "MFEMMaterial.h"


FEDomain* MFEMDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel& fem = *pmat->GetFEModel();
	FE_Element_Class eclass = spec.eclass;
	FEDomain* dom = nullptr;
   // MFEM Solver materials
   MFEMMaterial* mfem_mat = dynamic_cast<MFEMMaterial*>(pmat);
   MFEMUncoupledFEBioMaterial* uncoupled_mat = dynamic_cast<MFEMUncoupledFEBioMaterial*>(pmat);

   // FEBio Solver materials
   FEElasticMaterial* febio_elast_mat = dynamic_cast<FEElasticMaterial*>(pmat);
   FEUncoupledMaterial* febio_uncoupled_mat = dynamic_cast<FEUncoupledMaterial*>(pmat);

	if (mfem_mat || uncoupled_mat || febio_elast_mat || febio_uncoupled_mat)
	{
		if (eclass == FE_ELEM_SHELL)
		{
			dom = fecore_new<FEShellDomainOld>("mfem domain 2d", &fem);
		}
		else if (eclass == FE_ELEM_SOLID)
		{
         if (uncoupled_mat || mfem_mat)
         {
            // MFEM Solver
            dom = fecore_new<FESolidDomain>("mfem domain 3d", &fem);
         }
         else
         {
            // FEBio Solver
            if (febio_uncoupled_mat)
            {
               dom = fecore_new<FESolidDomain>("mfem immersed 3field febio domain", &fem);
            }
            else
            {
               dom = fecore_new<FESolidDomain>("mfem immersed febio domain", &fem);
            }            
         }
		}
		if (dom) dom->SetMaterial(pmat);
	}
	return dom;
}
