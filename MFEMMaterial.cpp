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