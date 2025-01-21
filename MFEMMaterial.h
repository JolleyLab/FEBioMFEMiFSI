#pragma once
#include <FECore/FEMaterial.h>

class MFEMMaterial : public FEMaterial
{
public:
	MFEMMaterial(FEModel* fem);

   // return material properties
   virtual double Density() = 0;
};

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
   double m_rho; // density
   double m_mu; // dynamic viscosity
   double m_kinvis; // kinematic viscosity (computed from density and dynamic viscosity)

	DECLARE_FECORE_CLASS();
};