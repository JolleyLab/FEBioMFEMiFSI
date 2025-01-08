#include <FECore/FECoreTask.h>

// Initial attempt to couple FEBio with the MFEM FSI solver using a task class. 
// The task will simply contruct an mfem mesh from the FEBio mesh and then 
// call the fluid solver.
class FEBioMFEMiFSITask : public FECoreTask
{
public:
	FEBioMFEMiFSITask(FEModel* fem) : FECoreTask(fem) {}

	// The Init function can be used to read a configuration file if necessary.
	// Currently, this is not used.
	bool Init(const char* szfile) override;

	// This will be called by FEBio after processing the input file. 
	bool Run() override;
};
