#pragma once
#include <FECore/FEInitialCondition.h>

class MFEMSolver;

// Base class for initial conditions that can be used in the MFEMSolver
class MFEMInitialCondition : public FEInitialCondition
{
public:
	MFEMInitialCondition(FEModel* fem) : FEInitialCondition(fem) {}

   virtual void CopyFrom(FEInitialCondition& ic)
   {
      MFEMInitialCondition* mic = dynamic_cast<MFEMInitialCondition*>(&ic);
      assert(mic);
      CopyParameterListState(mic->GetParameterList());
   };

	virtual void Apply(MFEMSolver* solver) = 0;
};

// Class that applies an initial pressure
class MFEMInitialPressure : public MFEMInitialCondition
{
public:
	MFEMInitialPressure(FEModel* fem) : MFEMInitialCondition(fem) {}

	void Apply(MFEMSolver* solver) override;

   void CopyFrom(FEInitialCondition& ic) override;

   // set the pressure
   void SetPressure(double p) { m_pressure = p; }

   // get the pressure
   double GetPressure() const;

private:
	FEParamDouble m_pressure;

	DECLARE_FECORE_CLASS();
};

// Class that applies an initial velocity
class MFEMInitialVelocity : public MFEMInitialCondition
{
public:
	MFEMInitialVelocity(FEModel* fem) : MFEMInitialCondition(fem) {}

	void Apply(MFEMSolver* solver) override;

   void CopyFrom(FEInitialCondition& ic) override;

private:
	FEParamDouble m_vx;
	FEParamDouble m_vy;
	FEParamDouble m_vz;

	DECLARE_FECORE_CLASS();
};
