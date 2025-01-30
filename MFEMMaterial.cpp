#include "MFEMMaterial.h"
#include "FECore/log.h"

MFEMMaterial::MFEMMaterial(FEModel* fem) : FEMaterial(fem)
{

}

BEGIN_FECORE_CLASS(MFEMNewtonianFluid, FEMaterial)
   ADD_PARAMETER(m_rho, FE_RANGE_GREATER(0.0), "rho")->setUnits(UNIT_DENSITY)->setLongName("density");
   ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu")->setUnits(UNIT_VISCOSITY)->setLongName("dynamic viscosity");
END_FECORE_CLASS();

MFEMNewtonianFluid::MFEMNewtonianFluid(FEModel* fem) : MFEMMaterial(fem)
{
   // defaults
   m_rho = 1.0;
   m_mu = 0.0;
   m_kinvis = 0.0;
}

double MFEMNewtonianFluid::Density()
{
   return m_rho;
}

double MFEMNewtonianFluid::DynamicViscosity()
{
   return m_mu;
}

double MFEMNewtonianFluid::KinematicViscosity()
{
   if (abs(m_rho) <= 0.0)
   {
      feLogError("Density is zero in MFEMNewtonianFluid");
   }
   return m_mu/m_rho;
}

/////////////////////////////////////////////////////////////////////////
// solid materials

MFEMSolidMaterial::MFEMSolidMaterial(FEModel* fem) : MFEMMaterial(fem)
{

}

BEGIN_FECORE_CLASS(MFEMCompressibleNeoHookeanMaterial, MFEMSolidMaterial)
   ADD_PARAMETER(m_rho, FE_RANGE_GREATER(0.0), "rho")->setUnits(UNIT_DENSITY)->setLongName("density");
   ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
   ADD_PARAMETER(m_lambda, FE_RANGE_GREATER_OR_EQUAL(0.0), "lambda");
END_FECORE_CLASS();

MFEMCompressibleNeoHookeanMaterial::MFEMCompressibleNeoHookeanMaterial(FEModel* fem) : MFEMSolidMaterial(fem)
{
   // defaults
   m_rho = 1.0;
   m_mu = 0.0;
   m_lambda = 0.0;

   m_rho_coeff = new FEScalarCoefficient(fem, m_rho);
   m_mu_coeff = new FEScalarCoefficient(fem, m_mu);
   m_lambda_coeff = new FEScalarCoefficient(fem, m_lambda);
}

MFEMNearlyIncompressibleMaterial::MFEMNearlyIncompressibleMaterial(FEModel* fem) : MFEMSolidMaterial(fem)
{

}

BEGIN_FECORE_CLASS(MFEMNearlyIncompressibleNeoHookeanMaterial, MFEMNearlyIncompressibleMaterial)
   ADD_PARAMETER(m_rho, FE_RANGE_GREATER(0.0), "rho")->setUnits(UNIT_DENSITY)->setLongName("density");
   ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
   ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(0.0), "kappa");
END_FECORE_CLASS();

MFEMNearlyIncompressibleNeoHookeanMaterial::MFEMNearlyIncompressibleNeoHookeanMaterial(FEModel* fem) : MFEMNearlyIncompressibleMaterial(fem)
{
   // defaults
   m_rho = 1.0;
   m_mu = 0.0;
   m_kappa = 0.0;

   m_rho_coeff = new FEScalarCoefficient(fem, m_rho);
   m_mu_coeff = new FEScalarCoefficient(fem, m_mu);
   m_kappa_coeff = new FEScalarCoefficient(fem, m_kappa);
}

BEGIN_FECORE_CLASS(MFEMNearlyIncompressibleAnnularDisk, MFEMNearlyIncompressibleMaterial)
   ADD_PARAMETER(m_rho, FE_RANGE_GREATER(0.0), "rho")->setUnits(UNIT_DENSITY)->setLongName("density");
   ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
   ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(0.0), "kappa");
   ADD_PARAMETER(m_w, "thickness")->setUnits(UNIT_LENGTH)->setLongName("disk radial thickness");
END_FECORE_CLASS();

MFEMNearlyIncompressibleAnnularDisk::MFEMNearlyIncompressibleAnnularDisk(FEModel* fem) : MFEMNearlyIncompressibleMaterial(fem)
{
   // defaults
   m_rho = 1.0;
   m_mu = 0.0;
   m_kappa = 0.0;
   m_w = 0.0;

   m_rho_coeff = new FEScalarCoefficient(fem, m_rho);
   m_mu_coeff = new FEScalarCoefficient(fem, m_mu);
   m_kappa_coeff = new FEScalarCoefficient(fem, m_kappa);
}