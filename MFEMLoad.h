#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelParam.h>

class MFEMSolver;

class MFEMLoad : public FESurfaceLoad
{
public:
	MFEMLoad(FEModel* fem) : FESurfaceLoad(fem) {}

	virtual void Apply(MFEMSolver* solver, int attr) = 0;
};

class MFEMNeumannBC : public MFEMLoad
{
public:
	MFEMNeumannBC(FEModel* fem) : MFEMLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:
	FEParamDouble m_tx;
	FEParamDouble m_ty;
	FEParamDouble m_tz;

	DECLARE_FECORE_CLASS();
};
