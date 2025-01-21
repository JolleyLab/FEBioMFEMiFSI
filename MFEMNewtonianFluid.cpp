#include "MFEMNewtonianFluid.h"
#include "FECore/log.h"

// only register the bottom level derived class, second arg is top level parent class 
BEGIN_FECORE_CLASS(MFEMNewtonianFluid, FEMaterial)
   ADD_PARAMETER(rho, FE_RANGE_GREATER(0.0), "rho")->setUnits(UNIT_DENSITY)->setLongName("density");
   ADD_PARAMETER(mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu")->setUnits(UNIT_VISCOSITY)->setLongName("dynamic viscosity");
END_FECORE_CLASS();

MFEMNewtonianFluid::MFEMNewtonianFluid(FEModel* fem) : MFEMMaterial(fem)
{
   rho = 0.0; // these will be default values in GUI
   mu = 0.0;
   kinvis = 0.0;
}

double MFEMNewtonianFluid::Density()
{
   return rho;
}

double MFEMNewtonianFluid::DynamicViscosity()
{
   return mu;
}

double MFEMNewtonianFluid::KinematicViscosity()
{
   if (abs(rho) <= 0.0)
   {
      feLogError("Density is zero in MFEMNewtonianFluid");
   }
   return mu/rho;
}