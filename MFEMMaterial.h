#pragma once
#include <FECore/FEMaterial.h>
#include "MFEMTools.h"
#include <FECore/log.h>
#include <Elasticity/Materials/nearly_incompressible_material.hpp>
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FECore/tens4d.h>
#include <FEBioMech/FEIncompNeoHookean.h>

class MFEMMaterial : public FEMaterial
{
public:
	MFEMMaterial(FEModel* fem);

   FECORE_BASE_CLASS(MFEMMaterial)
};

//-------------------------------------------------------------------------
/// MFEM Fluid Materials
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

//-------------------------------------------------------------------------
/// MFEM Solid Materials
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
   ~MFEMCompressibleNeoHookeanMaterial();

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
   ~MFEMNearlyIncompressibleNeoHookeanMaterial();

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
   ~MFEMNearlyIncompressibleAnnularDisk();

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

// Uncoupled materials using FEBio material classes
class MFEMUncoupledFEBioMaterial : public FEUncoupledMaterial, public NearlyIncompressibleMaterial
{
public:
   MFEMUncoupledFEBioMaterial(FEModel* fem);
   virtual ~MFEMUncoupledFEBioMaterial() noexcept;

   // set material point data
   void SetNearlyIncStrainMeasures(const mfem::Vector &u, const mfem::DenseMatrix &gradu) override;

   // Evaluate stress
   void EvalCauchyStress(mfem::DenseMatrix &stress) override;

   // Evaluate tangent
   double EvalElasticityTensor(int i, int j, int k, int l) override;

   // Set material constants, here we use it to map current gauss pt to physical space
   // for evaluating material parameters using FEParamDouble
   void SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T) override;

   FEScalarCoefficient *GetDensityCoeff() { return m_rho_coeff; }

private:
   FEMaterialPoint* m_materialpoint = nullptr;
   bool compute_tangent = true;
   tens4ds m_C;
   mfem::Vector r0;

   FEScalarCoefficient *m_rho_coeff = nullptr;

   FECORE_BASE_CLASS(MFEMUncoupledFEBioMaterial)
};

class FEBioMooneyRivlin : public MFEMUncoupledFEBioMaterial
{
public:
   FEBioMooneyRivlin(FEModel* pfem);
   ~FEBioMooneyRivlin() override {};

public:
	FEParamDouble	m_c1;	//!< Mooney-Rivlin coefficient C1
	FEParamDouble	m_c2;	//!< Mooney-Rivlin coefficient C2

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;   
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

class FEBioIncompNeoHookean : public MFEMUncoupledFEBioMaterial
{
public:
	FEBioIncompNeoHookean(FEModel* pfem);

public:
	FEParamDouble	m_G;	//!< Shear modulus

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// FEBio Materials

// FEBio version of MFEMNearlyIncompressibleAnnularDisk // special material for Annular Disk test problem
class FEAnnularDisk : public FEElasticMaterial
{
public:
   FEAnnularDisk(FEModel* pfem);

public:
   FEParamDouble m_mu; // shear modulus (c0) // zero for annular disk
   FEParamDouble m_kappa; // bulk modulus
   FEParamDouble m_w; // radial thickness

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;

   tens4ds prod_ikjl(const mat3ds&a, const mat3ds& b);
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};