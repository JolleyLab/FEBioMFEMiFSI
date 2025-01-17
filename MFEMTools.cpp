#include "MFEMTools.h"
#include <FECore/FEModel.h>

double FEScalarCoefficient::Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip)
{
	double coord[3];
	mfem::Vector xval(coord, 3);
	T.Transform(ip, xval);

	// the current time in mfem may not always match the current time in FEBio
	// so we have to do this little trick to make sure parameters are evaluated
	// at the correct time
	FETimeInfo& ti = m_fem->GetTime();
	double oldTime = ti.currentTime;
	ti.currentTime = GetTime();

	// setup the material point
	FEMaterialPoint pt;
	pt.m_r0 = vec3d(coord[0], coord[1], coord[2]);
	double v = m_val(pt);

	// reset the FEBio time
	ti.currentTime = oldTime;

	return v;
}

void FEVectorCoefficient::Eval(mfem::Vector& V, mfem::ElementTransformation& T,
	const mfem::IntegrationPoint& ip)
{
	double coord[3];
	mfem::Vector xval(coord, 3);
	T.Transform(ip, xval);

	// the current time in mfem may not always match the current time in FEBio
	// so we have to do this little trick to make sure parameters are evaluated
	// at the correct time
	FETimeInfo& ti = m_fem->GetTime();
	double oldTime = ti.currentTime;
	ti.currentTime = GetTime();

	// setup the material point
	FEMaterialPoint pt;
	pt.m_r0 = vec3d(coord[0], coord[1], coord[2]);
	double vx = m_vx(pt);
	double vy = m_vy(pt);
	double vz = m_vz(pt);

	// reset the FEBio time
	ti.currentTime = oldTime;

	V.SetSize(3);
	V(0) = vx;
	V(1) = vy;
	V(2) = vz;
}
