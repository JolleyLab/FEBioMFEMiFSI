#include "MFEMCallbacks.h"

#include <FECore/FEAnalysis.h>
#include <FECore/FESolver.h>

#include "../MFEMCore/MFEMSolver.h"

bool timing_callback(FEModel* fem, unsigned int nreason, void* pd)
{
   FEAnalysis* analysis = fem->GetCurrentStep();
   FESolver* solver = analysis->GetFESolver();
   auto mfem_solver = dynamic_cast<MFEMSolver*>(solver);
   MFEM_ASSERT(mfem_solver, "Error: not a MFEM solver");
   mfem_solver->PostSolve();
   return true;
}

bool time_step_callback(FEModel* fem, unsigned int nreason, void* pd)
{
   FEAnalysis* analysis = fem->GetCurrentStep();
   FESolver* solver = analysis->GetFESolver();
   auto mfem_solver = dynamic_cast<MFEMSolver*>(solver);
   MFEM_ASSERT(mfem_solver, "Error: not a MFEM solver");
   mfem_solver->AdjustTimeStep();
   return true;
}