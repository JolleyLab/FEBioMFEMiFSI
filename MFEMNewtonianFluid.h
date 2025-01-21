#pragma once
#include "MFEMMaterial.h"

class MFEMNewtonianFluid : public MFEMMaterial
{
public:
	MFEMNewtonianFluid(FEModel* fem);

   // return material properties
   double Density() override;
   double DynamicViscosity();
   double KinematicViscosity();

public:
   // TO DO: Think about making these variable, but for now assume constant
   double rho; // density
   double mu; // dynamic viscosity
   double kinvis; // kinematic viscosity (computed from density and dynamic viscosity)

	DECLARE_FECORE_CLASS();
};
