#include "MFEMDomain.h"
#include <FECore/FEModel.h>
#include "MFEMMaterial.h"
#include <FECore/FEAnalysis.h>
#include "MFEMSolver.h"

FEDomain* MFEMDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel& fem = *pmat->GetFEModel();
   FE_Element_Class eclass = spec.eclass;
	FEDomain* dom = nullptr;

   FESolver &solver = *fem.GetCurrentStep()->GetFESolver();
   MFEMSolver* mfem_solver = dynamic_cast<MFEMSolver*>(&solver);
   if (mfem_solver == nullptr)
   {
      MFEM_ABORT("Invalid solver type");
   }
   bool immersed_solid = mfem_solver->ImmersedSolid();
   int solidsolvertype = mfem_solver->SolidSolverType();	

   MFEMMaterial* mfem_mat = dynamic_cast<MFEMMaterial*>(pmat); // MFEM materials
   FEElasticMaterial* febio_elast_mat = dynamic_cast<FEElasticMaterial*>(pmat); // FEBio materials
	if (mfem_mat || febio_elast_mat)
	{
		if (eclass == FE_ELEM_SHELL)
		{
			dom = fecore_new<FEShellDomainOld>("mfem domain 2d", &fem);
		}
		else if (eclass == FE_ELEM_SOLID)
		{
         if (mfem_mat)
         {
            // MFEM Material + MFEM Solver
            dom = fecore_new<FESolidDomain>("mfem domain 3d", &fem);
         }
         else if (febio_elast_mat)
         {
            if (solidsolvertype == MFEMSolver::SolidSolverType::FEBio)
            {
               if (dynamic_cast<FEUncoupledMaterial*>(pmat))
               {
                  // FEBio 3Field Domain + FEBio Material
                  dom = fecore_new<FESolidDomain>("mfem immersed 3field febio domain", &fem);
               }
               else
               {
                  // FEBio Domain + FEBio Material
                  dom = fecore_new<FESolidDomain>("mfem immersed febio domain", &fem);
               }   
            }
            else if (solidsolvertype == MFEMSolver::SolidSolverType::MFEM)
            {
               // MFEM Solver + FEBio Material
               dom = fecore_new<FESolidDomain>("mfem domain 3d", &fem);
            }
            else
            {
               MFEM_ABORT("Invalid solid solver type");
            }   
         }
         else
         {
            MFEM_ABORT("Invalid material");
         }
		}
      else
      {
         MFEM_ABORT("Unsupported element class");
      }
		if (dom) dom->SetMaterial(pmat);
	}
	return dom;
}
