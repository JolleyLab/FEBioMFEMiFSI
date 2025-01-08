#include <FECore/sdk.h>
#include <FECore/FEModule.h>
#include <FECore/FEModel.h>
#include "FEBioMFEMiFSITask.h"

class FEMFEMiFSI : public FEModule
{
public:
	FEMFEMiFSI() {}
	void InitModel(FEModel* fem)
	{
		DOFS& dofs = fem->GetDOFS();
		int pressure = dofs.AddVariable("pressure");
		dofs.SetDOFName(pressure, 0, "p");
		int velocity = dofs.AddVariable("velocity", VAR_VEC2);
		dofs.SetDOFName(velocity, 0, "vx");
		dofs.SetDOFName(velocity, 1, "vy");
		dofs.SetDOFName(velocity, 2, "vz");
	}
};

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

FECORE_PLUGIN void GetPluginVersion(int& major, int& minor, int& patch)
{
	major = 2;
	minor = 0;
	patch = 0;
}

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);

	// create the heat-transfer module
	const char* info = \
		"{ "
		"   \"title\" : \"MFM-iFSI\","
		"   \"info\"  : \"Interface to MFEM based iFSI solver\","
		"   \"author\": \"Ryan Black, Steve Maas\","
		"   \"version\": \"1.0\""
		"}";

	//	fecore.CreateModule(new FEMFEMiFSI, "mfemifsi", info);

	REGISTER_FECORE_CLASS(FEBioMFEMiFSITask, "mfemifsi");
}
