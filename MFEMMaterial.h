#pragma once
#include <FECore/FEMaterial.h>
#include "MFEMTools.h"
#include <FECore/log.h>

class MFEMMaterial : public FEMaterial
{
public:
	MFEMMaterial(FEModel* fem);
};

class MFEMNewtonianFluid : public MFEMMaterial
{
public:
	MFEMNewtonianFluid(FEModel* fem);

   // return material properties
   double Density();
   double DynamicViscosity();
   double KinematicViscosity();

private:
   // TO DO: Think about making these variable, but for now assume constant
   double m_rho; // density
   double m_mu; // dynamic viscosity
   double m_kinvis; // kinematic viscosity (computed from density and dynamic viscosity)

	DECLARE_FECORE_CLASS();
};

// solid materials
class MFEMSolidMaterial : public MFEMMaterial
{
public:
   MFEMSolidMaterial(FEModel* fem);

   virtual FEScalarCoefficient *GetDensityCoeff() = 0;
   virtual FEScalarCoefficient *GetC0Coeff() = 0;
};

class MFEMCompressibleNeoHookeanMaterial : public MFEMSolidMaterial
{
public:
   MFEMCompressibleNeoHookeanMaterial(FEModel* fem);

   FEScalarCoefficient *GetDensityCoeff() override { return m_rho_coeff; }
   FEScalarCoefficient *GetC0Coeff() override{ return m_mu_coeff; }
   FEScalarCoefficient *GetC1Coeff() { return m_lambda_coeff; }

private:
   FEParamDouble m_rho; // density
   FEParamDouble m_mu; // shear modulus (c0)
   FEParamDouble m_lambda; // lambda (c1)

   FEScalarCoefficient *m_rho_coeff = nullptr;
   FEScalarCoefficient *m_mu_coeff = nullptr;
   FEScalarCoefficient *m_lambda_coeff = nullptr;

   DECLARE_FECORE_CLASS();
};

class MFEMNearlyIncompressibleMaterial : public MFEMSolidMaterial
{
public:
   MFEMNearlyIncompressibleMaterial(FEModel* fem);

   virtual FEScalarCoefficient *GetBulkModulus() = 0;
};

class MFEMNearlyIncompressibleNeoHookeanMaterial : public MFEMNearlyIncompressibleMaterial
{
public:
   MFEMNearlyIncompressibleNeoHookeanMaterial(FEModel* fem);

   FEScalarCoefficient *GetDensityCoeff() override { return m_rho_coeff; }
   FEScalarCoefficient *GetC0Coeff() override{ return m_mu_coeff; }
   FEScalarCoefficient *GetBulkModulus() override { return m_kappa_coeff; }

private:
   FEParamDouble m_rho; // density
   FEParamDouble m_mu; // shear modulus (c0)
   FEParamDouble m_kappa; // bulk modulus

   FEScalarCoefficient *m_rho_coeff = nullptr;
   FEScalarCoefficient *m_mu_coeff = nullptr;
   FEScalarCoefficient *m_kappa_coeff = nullptr;

   DECLARE_FECORE_CLASS();
};

// special material for Annular Disk test problem
class MFEMNearlyIncompressibleAnnularDisk : public MFEMNearlyIncompressibleMaterial
{
public:
   MFEMNearlyIncompressibleAnnularDisk(FEModel* fem);

   FEScalarCoefficient *GetDensityCoeff() override { return m_rho_coeff; }
   FEScalarCoefficient *GetC0Coeff() override{ return m_mu_coeff; }
   FEScalarCoefficient *GetBulkModulus() override { return m_kappa_coeff; }
   double GetRadialThickness() 
   {
      if(m_w.isConst() == false)
      {
         feLogError("Radial thickness must be constant");
         assert(false);
      }
      return m_w.constValue(); 
   }

private:
   FEParamDouble m_rho; // density
   FEParamDouble m_mu; // shear modulus (c0)
   FEParamDouble m_kappa; // bulk modulus
   FEParamDouble m_w; // radial thickness

   FEScalarCoefficient *m_rho_coeff = nullptr;
   FEScalarCoefficient *m_mu_coeff = nullptr;
   FEScalarCoefficient *m_kappa_coeff = nullptr;

   DECLARE_FECORE_CLASS();
};