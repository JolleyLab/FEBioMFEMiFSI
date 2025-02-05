#include "MFEMMaterial.h"
#include "FECore/log.h"
#include <FECore/tens4d.h>
#include <FECore/mat3d.h>
#include <Elasticity/Materials/material_type.hpp>
#include <FECore/FEModel.h>

MFEMMaterial::MFEMMaterial(FEModel* fem) : FEMaterial(fem)
{

}

BEGIN_FECORE_CLASS(MFEMNewtonianFluid, MFEMMaterial)
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

MFEMCompressibleNeoHookeanMaterial::~MFEMCompressibleNeoHookeanMaterial()
{
   delete m_rho_coeff;
   delete m_mu_coeff;
   delete m_lambda_coeff;
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

MFEMNearlyIncompressibleNeoHookeanMaterial::~MFEMNearlyIncompressibleNeoHookeanMaterial()
{
   delete m_rho_coeff;
   delete m_mu_coeff;
   delete m_kappa_coeff;
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

MFEMNearlyIncompressibleAnnularDisk::~MFEMNearlyIncompressibleAnnularDisk()
{
   delete m_rho_coeff;
   delete m_mu_coeff;
   delete m_kappa_coeff;
}

MFEMUncoupledFEBioMaterial::MFEMUncoupledFEBioMaterial(FEModel* fem) : FEUncoupledMaterial(fem), NearlyIncompressibleMaterial(mat::MaterialType::UpdatedLagrangian,3)
{
   m_rho_coeff = new FEScalarCoefficient(fem, m_density);

   // init to zero
   m_C.zero();
}

MFEMUncoupledFEBioMaterial::~MFEMUncoupledFEBioMaterial() noexcept
{
   delete m_materialpoint;
   delete m_rho_coeff;
}

void MFEMUncoupledFEBioMaterial::SetNearlyIncStrainMeasures(const mfem::Vector &u, const mfem::DenseMatrix &gradu)
{
   // setup material point
   if(m_materialpoint) delete m_materialpoint;
   m_materialpoint = new FEMaterialPoint((*this).CreateMaterialPointData()); // this will delete the ElasticMaterialPoint data (m_data) as well

   // set relevant data
   FEElasticMaterialPoint& pt = *m_materialpoint->ExtractData<FEElasticMaterialPoint>(); // this references m_data (MaterialPointData) in m_materialpoint (MaterialPoint)
   pt.m_J = J;
   
   pt.m_F(0,0) = F(0,0);   
   pt.m_F(0,1) = F(0,1);   
   pt.m_F(0,2) = F(0,2);
   pt.m_F(1,0) = F(1,0);
   pt.m_F(1,1) = F(1,1);
   pt.m_F(1,2) = F(1,2);
   pt.m_F(2,0) = F(2,0);
   pt.m_F(2,1) = F(2,1);
   pt.m_F(2,2) = F(2,2);

   pt.m_buncoupled = true;

   m_materialpoint->m_r0 = vec3d(r0(0),r0(1),r0(2));

   // flag to recompute tangent at new gauss pt
   compute_tangent = true;
}

void MFEMUncoupledFEBioMaterial::SetMaterialConstants(const mfem::IntegrationPoint &ip, mfem::ElementTransformation &T)
{
   // note this function is called prior to SetNearlyIncStrainMeasures() above
   // here we use it to get the coordinates of the current integration point
   r0.SetSize(3);
	T.Transform(ip, r0);
}

void MFEMUncoupledFEBioMaterial::EvalCauchyStress(DenseMatrix &stress)
{
   stress.SetSize(vdim);

   // Compute sigma = sigma_iso + sigma_vol
   mat3ds m_s = (*this).Stress(*m_materialpoint);

   // copy data to mfem stress matrix
   stress(0,0) = m_s(0,0);
   stress(0,1) = m_s(0,1);
   stress(0,2) = m_s(0,2);
   stress(1,0) = m_s(1,0);
   stress(1,1) = m_s(1,1);
   stress(1,2) = m_s(1,2);
   stress(2,0) = m_s(2,0);
   stress(2,1) = m_s(2,1);
   stress(2,2) = m_s(2,2);
}

double MFEMUncoupledFEBioMaterial::EvalElasticityTensor(int i, int j, int k, int l)
{
   double val = 0.0;
   if (formulation_type == form_type::UpdatedLagrangian)
   {
      if (compute_tangent) 
      {
         m_C = (*this).Tangent(*m_materialpoint);
         compute_tangent = false;
      }
      val = m_C(i,j,k,l);
   }
   else
   {
      MFEM_ABORT("Error: TotalLagrangian formulation not implemented");
   }
   
   return val;
}

BEGIN_FECORE_CLASS(FEBioMooneyRivlin, MFEMUncoupledFEBioMaterial)
   ADD_PARAMETER(m_c1, "c1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c2, "c2")->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

FEBioMooneyRivlin::FEBioMooneyRivlin(FEModel* pfem) : MFEMUncoupledFEBioMaterial(pfem)
{
   m_c1 = 0.0;
   m_c2 = 0.0;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric stress
mat3ds FEBioMooneyRivlin::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get material parameters
	double c1 = m_c1(mp);
	double c2 = m_c2(mp);

	// determinant of deformation gradient
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	//
	// W = C1*(I1 - 3) + C2*(I2 - 3)
	//
	// Wi = dW/dIi
	double W1 = c1;
	double W2 = c2;
	// ---

	// calculate T = F*dW/dC*Ft
	// T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric tangent
tens4ds FEBioMooneyRivlin::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get material parameters
	double c1 = m_c1(mp);
	double c2 = m_c2(mp);

	// determinant of deformation gradient
	double J = pt.m_J;
	double Ji = 1.0/J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1, W2;
	W1 = c1;
	W2 = c2;
	// ---

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds T = B*(W1 + W2*I1) - B2*W2;
	mat3ds devs = T.dev()*(2.0/J);

	// Identity tensor
	mat3ds I(1,1,1,0,0,0);

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// d2W/dCdC:C
	mat3ds WCCxC = B*(W2*I1) - B2*W2;

	tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);

	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c;
}

// define the material parameters
BEGIN_FECORE_CLASS(FEBioIncompNeoHookean, MFEMUncoupledFEBioMaterial)
	ADD_PARAMETER(m_G, FE_RANGE_GREATER(0.0), "G")->setUnits(UNIT_PRESSURE)->setLongName("shear modulus");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEBioIncompNeoHookean::FEBioIncompNeoHookean(FEModel* pfem) : MFEMUncoupledFEBioMaterial(pfem)
{
	m_G = 0.0;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric stress
mat3ds FEBioIncompNeoHookean::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// shear modulus
	double G = m_G(mp);

	// deformation gradient
	double J = pt.m_J;

	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();

	// calculate deviatoric stress
	return (B - mat3dd(B.tr()/3.))*(G/pow(J, 5.0/3.0));
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FEBioIncompNeoHookean::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// shear modulus
	double G = m_G(mp);

	// deformation gradient
	double J = pt.m_J;

	// left cauchy-green matrix (i.e. the 'b' matrix)
	mat3ds B = pt.LeftCauchyGreen();

	// trace of b
	double Ib = B.tr();

	double muJ = G/pow(J, 5.0/3.0);

	mat3ds I(1,1,1,0,0,0);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxI = dyad1s(B, I); // = BxI + IxB

	return (I4*Ib -BxI + IxI*(Ib/3))*(2.0*muJ/3.0);
}