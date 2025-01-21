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

class MFEMTractionFreeBC : public MFEMLoad
{
public:
	MFEMTractionFreeBC(FEModel* fem) : MFEMLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:

	DECLARE_FECORE_CLASS();
};

class MFEMPressureBC : public MFEMLoad
{
public:
	MFEMPressureBC(FEModel* fem) : MFEMLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:
   // applies p*nor to the surface
	FEParamDouble m_p;

	DECLARE_FECORE_CLASS();
};

class MFEMBackflowStabilizationBC : public MFEMLoad
{
public:
	MFEMBackflowStabilizationBC(FEModel* fem);

	void Apply(MFEMSolver* solver, int attr) override;

private:
   // assumes constant beta
	FEParamDouble m_beta;

	DECLARE_FECORE_CLASS();
};

class MFEMResistanceBC : public MFEMLoad
{
public:
	MFEMResistanceBC(FEModel* fem) : MFEMLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:
   // assumes constant R and p0
	FEParamDouble m_R; // proximal resistance
	FEParamDouble m_p0; // initial pressure

	DECLARE_FECORE_CLASS();
};

class MFEMWindkesselBC : public MFEMLoad
{
public:
	MFEMWindkesselBC(FEModel* fem) : MFEMLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

private:
   // assumes constant R, Rd, C, p0(spatial), pd
	FEParamDouble m_R; // proximal resistance 
	FEParamDouble m_Rd; // distal resistance
   FEParamDouble m_C; // vascular compliance
   FEParamDouble m_p0; // initial pressure // TODO make this spatially dependent thru a coeff
   FEParamDouble m_pd; // distal pressure
   FEParamDouble m_Q0; // initial flow rate

	DECLARE_FECORE_CLASS();
};

