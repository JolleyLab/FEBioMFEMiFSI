#include "MFEMTools.h"
#include <FECore/FEModel.h>

double FEScalarCoefficient::Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip)
{
	double coord[3];
	mfem::Vector xval(coord, 3);
	T.Transform(ip, xval);

	// setup the material point
	FEMaterialPoint pt;
	pt.m_r0 = vec3d(coord[0], coord[1], coord[2]);
	double v = m_val(pt);

	return v;
}

void FEVectorCoefficient::Eval(mfem::Vector& V, mfem::ElementTransformation& T,
	const mfem::IntegrationPoint& ip)
{
	double coord[3];
	mfem::Vector xval(coord, 3);
	T.Transform(ip, xval);

	// setup the material point
	FEMaterialPoint pt;
	pt.m_r0 = vec3d(coord[0], coord[1], coord[2]);
	double vx = m_vx(pt);
	double vy = m_vy(pt);
	double vz = m_vz(pt);

	V.SetSize(3);
	V(0) = vx;
	V(1) = vy;
	V(2) = vz;
}