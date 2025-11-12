#include "MFEMMaterial.h"
#include "FECore/log.h"
#include <FECore/tens4d.h>
#include <FECore/mat3d.h>
#include "../MFEMFSI/Elasticity/Materials/material_type.hpp"
#include <FECore/FEModel.h>

// using namespace mfem;

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

BEGIN_FECORE_CLASS(FEAnnularDisk, FEElasticMaterial)
   ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
   ADD_PARAMETER(m_w, "thickness")->setUnits(UNIT_LENGTH)->setLongName("disk radial thickness");
END_FECORE_CLASS();

FEAnnularDisk::FEAnnularDisk(FEModel* pfem) : FEElasticMaterial(pfem) 
{
	m_mu = 0.0;
   m_w = 0.0;
   m_kappa = 0.0;
}

//-----------------------------------------------------------------------------
mat3ds FEAnnularDisk::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

   // compute ehat_theta
   vec3d curr_pos = mp.m_rt;
   // shift to center of disk
   curr_pos(0) -= 0.5;
   curr_pos(1) -= 0.5;
   double mag = curr_pos.unit();
   vec3d ehat_theta(-1.0*curr_pos(1),curr_pos(0),0.0);   

   // shear modulus
   double mu_val = m_mu(mp);

   mat3d F = pt.m_F;
   vec3d Ftheta = F*ehat_theta;
   mat3d FtFt_temp = Ftheta & Ftheta; // this is sym
   mat3ds FtFt = FtFt_temp.sym(); // to change type

   return FtFt*(mu_val/pt.m_J);
}

//-----------------------------------------------------------------------------
tens4ds FEAnnularDisk::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

   // compute ehat_theta
   vec3d curr_pos = mp.m_rt;
   // shift to center of disk
   curr_pos(0) -= 0.5;
   curr_pos(1) -= 0.5;
   double mag = curr_pos.unit();
   vec3d ehat_theta(-1.0*curr_pos(1),curr_pos(0),0.0);   

   // shear modulus
   double mu_val = m_mu(mp);

   mat3d F = pt.m_F;
   vec3d Ftheta = F*ehat_theta;
   mat3d FtFt_temp = Ftheta & Ftheta; // this is sym
   mat3ds FtFt = FtFt_temp.sym(); // to change type

   // left cauchy-green matrix (i.e. the 'b' matrix)
	mat3ds b = pt.LeftCauchyGreen();

   // c_ijkl = mu/J * (b_ik FtFt_jl)
   tens4ds c = prod_ikjl(b, FtFt);
	return c*(mu_val/pt.m_J);
}

//()_ijkl = a_ik * b_jl
tens4ds FEAnnularDisk::prod_ikjl(const mat3ds&a, const mat3ds& b)
{
	tens4ds c;

    c.d[ 0] = a.xx() * b.xx(); // 0000
    
    c.d[ 1] = a.xy() * b.xy(); // 0011
    c.d[ 2] = a.yy() * b.yy(); // 1111

    c.d[ 3] = a.xz() * b.xz(); // 0022
    c.d[ 4] = a.yz() * b.yz(); // 1122
    c.d[ 5] = a.zz() * b.zz(); // 2222
    
    c.d[ 6] = a.xx() * b.xy(); // 0001
    c.d[ 7] = a.xy() * b.yy(); // 1101
    c.d[ 8] = a.xz() * b.yz(); // 2201
    c.d[ 9] = a.xx() * b.yy(); // 0101

    c.d[10] = a.xy() * b.xz(); // 0012
    c.d[11] = a.yy() * b.yz(); // 1112
    c.d[12] = a.yz() * b.zz(); // 2212
    c.d[13] = a.xy() * b.yz(); // 0112
    c.d[14] = a.yy() * b.zz(); // 1212

    c.d[15] = a.xx() * b.xz(); // 0002
    c.d[16] = a.xy() * b.yz(); // 1102
    c.d[17] = a.xz() * b.zz(); // 2202
    c.d[18] = a.xx() * b.yz(); // 0102
    c.d[19] = a.xy() * b.zz(); // 1202
    c.d[20] = a.xx() * b.zz(); // 0202

	return c;
}


//-----------------------------------------------------------------------------
double FEAnnularDisk::StrainEnergyDensity(FEMaterialPoint& mp)
{
   MFEM_ABORT("Error: StrainEnergyDensity not implemented");
	return 0.0;
}

//-----------------------------------------------------------------------------
FEBioElasticMaterial::FEBioElasticMaterial(FEElasticMaterial* mat_) : emat(mat_), MaterialType(mat::MaterialType::UpdatedLagrangian, 3)
{
	FEParam* p = emat->GetParameter("density"); assert(p);
	FEParamDouble& pd = p->value<FEParamDouble>();
	m_rho_coeff = new FEScalarCoefficient(emat->GetFEModel(), pd);
}

void FEBioElasticMaterial::SetMaterialConstants(const mfem::IntegrationPoint& ip, mfem::ElementTransformation& T)
{
	// note this function is called prior to SetAdditionalData()
	// here we use it to get the coordinates of the current integration point
	r0.SetSize(3);
	T.Transform(ip, r0);
}

void FEBioElasticMaterial::SetAdditionalData(const mfem::Vector& u, const mfem::DenseMatrix& gradu)
{
	// setup material point
	if (m_materialpoint) delete m_materialpoint;
	m_materialpoint = new FEMaterialPoint(emat->CreateMaterialPointData()); // this will delete the ElasticMaterialPoint data (m_data) as well

	// set relevant data
	FEElasticMaterialPoint& pt = *m_materialpoint->ExtractData<FEElasticMaterialPoint>(); // this references m_data (MaterialPointData) in m_materialpoint (MaterialPoint)
	pt.m_J = J;

	pt.m_F(0, 0) = F(0, 0);
	pt.m_F(0, 1) = F(0, 1);
	pt.m_F(0, 2) = F(0, 2);
	pt.m_F(1, 0) = F(1, 0);
	pt.m_F(1, 1) = F(1, 1);
	pt.m_F(1, 2) = F(1, 2);
	pt.m_F(2, 0) = F(2, 0);
	pt.m_F(2, 1) = F(2, 1);
	pt.m_F(2, 2) = F(2, 2);

	m_materialpoint->m_r0 = vec3d(r0(0), r0(1), r0(2));

	// flag to recompute tangent at new gauss pt
	compute_tangent = true;
}

// Evaluate stress
void FEBioElasticMaterial::EvalCauchyStress(mfem::DenseMatrix& stress)
{
	stress.SetSize(vdim);

	// Compute sigma
	mat3ds m_s = emat->Stress(*m_materialpoint);

	// copy data to mfem stress matrix
	stress(0, 0) = m_s(0, 0);
	stress(0, 1) = m_s(0, 1);
	stress(0, 2) = m_s(0, 2);
	stress(1, 0) = m_s(1, 0);
	stress(1, 1) = m_s(1, 1);
	stress(1, 2) = m_s(1, 2);
	stress(2, 0) = m_s(2, 0);
	stress(2, 1) = m_s(2, 1);
	stress(2, 2) = m_s(2, 2);
}

// Evaluate tangent
double FEBioElasticMaterial::EvalElasticityTensor(int i, int j, int k, int l)
{
	double val = 0.0;
	if (formulation_type == form_type::UpdatedLagrangian)
	{
		if (compute_tangent)
		{
			m_C = emat->Tangent(*m_materialpoint);
			compute_tangent = false;
		}
		val = m_C(i, j, k, l);
	}
	else
	{
		MFEM_ABORT("Error: TotalLagrangian formulation not implemented");
	}

	return val;
}
