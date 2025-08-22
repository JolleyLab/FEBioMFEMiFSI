#include <FECore/sdk.h>
#include <FECore/FEModule.h>
#include <FECore/FEModel.h>
#include <FEBioMech/FEIncompNeoHookean.h>
#include <FEBioMech/FENeoHookean.h>
#include <FEBioMech/FEStVenantKirchhoff.h>
#include <FEBioMech/FEFungOrthotropic.h>

#include "MFEMAnalysis.h"
#include "MFEMDomain.h"
#include "MFEMSolver.h"
#include "MFEMMaterial.h"
#include "MFEMPlotVariables.h"
#include "MFEMInitialCondition.h"
#include "MFEMBoundaryCondition.h"
#include "MFEMLoad.h"
#include "MFEMImmersedElasticDomain.h"
#include "MFEMImmersed3FieldElasticDomain.h"
#include "mfem.hpp"

class FEMFEMiFSI : public FEModule
{
public:
	FEMFEMiFSI() {}
	void InitModel(FEModel* fem)
	{
		DOFS& dofs = fem->GetDOFS();
		int pressure = dofs.AddVariable("pressure");
		dofs.SetDOFName(pressure, 0, "p");
		int velocity = dofs.AddVariable("velocity", VAR_VEC3);
		dofs.SetDOFName(velocity, 0, "vx");
		dofs.SetDOFName(velocity, 1, "vy");
		dofs.SetDOFName(velocity, 2, "vz");

      // immersed solid
      int immersed_displacement = dofs.AddVariable("immersed_displacement", VAR_VEC3);
      dofs.SetDOFName(immersed_displacement, 0, "im_ux");
      dofs.SetDOFName(immersed_displacement, 1, "im_uy");
      dofs.SetDOFName(immersed_displacement, 2, "im_uz");
      int immersed_velocity = dofs.AddVariable("immersed_velocity", VAR_VEC3);
      dofs.SetDOFName(immersed_velocity, 0, "im_vx");
      dofs.SetDOFName(immersed_velocity, 1, "im_vy");
      dofs.SetDOFName(immersed_velocity, 2, "im_vz");
	}
};

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

FECORE_PLUGIN void GetPluginVersion(int& major, int& minor, int& patch)
{
	major = 1;
	minor = 0;
	patch = 0;
}

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
   FECoreKernel::SetInstance(&fecore);

   #ifdef MFEM_USE_MPI
   // initialize MPI and Hypre libraries
   
   // int argc;
   // char** argv;
   // fecore.GetCommandLineInput(argc,argv);
   // mfem::Mpi::Init(argc, argv);

   // works for most MPI implementations
   mfem::Mpi::Init(nullptr, nullptr);
   
   mfem::Hypre::Init();
   #endif	

	// create the new module
	const char* info = \
		"{ "
		"   \"title\" : \"MFEM-iFSI\","
		"   \"info\"  : \"Interface to MFEM based iFSI solver\","
		"   \"author\": \"Ryan Black, Steve Maas\","
		"   \"version\": \"1.0\""
		"}";

	fecore.CreateModule(new FEMFEMiFSI, "mfemifsi", info);
   // fecore.AddModuleDependency("solid");

	// analyis classes (default type must match module name!)
	REGISTER_FECORE_CLASS(MFEMAnalysis, "mfemifsi");

	// Domain factory
	fecore.RegisterDomain(new MFEMDomainFactory);
	REGISTER_FECORE_CLASS(MFEMDomain2D, "mfem domain 2d");
	REGISTER_FECORE_CLASS(MFEMSolidDomain, "mfem domain 3d");
   REGISTER_FECORE_CLASS(MFEMImmersedElasticDomain, "mfem immersed febio domain");
   REGISTER_FECORE_CLASS(MFEMImmersed3FieldElasticDomain, "mfem immersed 3field febio domain");

	// feature classes
	REGISTER_FECORE_CLASS(MFEMSolver, "mfem");
   REGISTER_FECORE_CLASS(MFEMNewtonianFluid, "mfem newtonian fluid");

	REGISTER_FECORE_CLASS(MFEMInitialPressure, "initial pressure");
	REGISTER_FECORE_CLASS(MFEMInitialVelocity, "initial velocity");

   // velocity bcs
	REGISTER_FECORE_CLASS(MFEMVelDirichletBC, "prescribed velocity");
   REGISTER_FECORE_CLASS(MFEMNoSlipBC, "no slip");
   REGISTER_FECORE_CLASS(MFEMPrescribedVelocityComponentBC, "prescribed velocity component");

   // traction bcs
	REGISTER_FECORE_CLASS(MFEMNeumannBC, "traction");
   REGISTER_FECORE_CLASS(MFEMTractionFreeBC, "traction free");
   REGISTER_FECORE_CLASS(MFEMPressureBC, "pressure");

   // periodic bcs
   REGISTER_FECORE_CLASS(MFEMPeriodicBC, "periodic");
   
   // cardiovascular bcs
   // REGISTER_FECORE_CLASS(MFEMBackflowStabilizationBC, "backflow stabilization");
   REGISTER_FECORE_CLASS(MFEMResistanceBC, "resistance (R)");
   REGISTER_FECORE_CLASS(MFEMWindkesselBC, "windkessel (RCR)");
   REGISTER_FECORE_CLASS(MFEMMappedInflowBC, "mapped inflow");

   // body force
   REGISTER_FECORE_CLASS(MFEMBodyForce, "body force");

   // Set pressure nodal value
   REGISTER_FECORE_CLASS(MFEMSetNodalPressure, "set nodal pressure");

   // Immersed solid MFEM materials
   REGISTER_FECORE_CLASS(MFEMCompressibleNeoHookeanMaterial, "mfem compressible neo-hookean");
   REGISTER_FECORE_CLASS(MFEMNearlyIncompressibleNeoHookeanMaterial, "mfem nearly-incompressible neo-hookean");
   REGISTER_FECORE_CLASS(MFEMNearlyIncompressibleAnnularDisk, "mfem nearly-incompressible annular disk");

   // Immersed solid FEBio materials
   REGISTER_FECORE_CLASS(FEAnnularDisk, "annular disk");
   REGISTER_FECORE_CLASS(FEIncompNeoHookean, "incompressible neo-hookean");
   REGISTER_FECORE_CLASS(FENeoHookean, "neo-hookean");
   REGISTER_FECORE_CLASS(FEStVenantKirchhoff, "st. venant-kirchhoff");
   REGISTER_FECORE_CLASS(FEFungOrthotropic, "fung orthotropic");

   // Immersed solid bcs
   // Surface Loads
   REGISTER_FECORE_CLASS(ImmersedTractionLoad, "immersed traction");
   REGISTER_FECORE_CLASS(ImmersedPressureLoad, "immersed pressure load");
   REGISTER_FECORE_CLASS(ImmersedPresDispBC, "immersed prescribed displacement");

   // plot variables
	REGISTER_FECORE_CLASS(FEPltMFemPressure, "mfem pressure");
	REGISTER_FECORE_CLASS(FEPltMFemVelocity, "mfem velocity");
   REGISTER_FECORE_CLASS(FEPltMFemImmersedDisplacement, "mfem immersed displacement");
	REGISTER_FECORE_CLASS(FEPltMFemImmersedVelocity, "mfem immersed velocity");
}
