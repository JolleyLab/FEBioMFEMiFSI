#pragma once
#include <FECore/FEModelParam.h>
#include <FECore/FEMaterialPoint.h>
#include <mfem.hpp>

class FEScalarCoefficient : public mfem::Coefficient
{
private:
	FEModel* m_fem = nullptr;
	FEParamDouble& m_val;

public:
	FEScalarCoefficient(FEModel* fem, FEParamDouble& val) : m_fem(fem), m_val(val) {}

	double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override;
};

class FEVectorCoefficient : public mfem::VectorCoefficient
{
private:
	FEModel* m_fem = nullptr;
	FEParamDouble& m_vx;
	FEParamDouble& m_vy;
	FEParamDouble& m_vz;

public:
	FEVectorCoefficient(FEModel* fem, FEParamDouble& vx, FEParamDouble& vy, FEParamDouble& vz) : m_fem(fem), m_vx(vx), m_vy(vy), m_vz(vz), VectorCoefficient(3) {}

	void Eval(mfem::Vector& V, mfem::ElementTransformation& T,
		const mfem::IntegrationPoint& ip) override;
};
