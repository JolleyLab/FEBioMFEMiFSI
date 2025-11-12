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

// Class that applies initial conditions from a folder
// written in MFEM mesh format (see option write_final_state_as_IC)
class MFEMInitialConditionFromFolder : public MFEMInitialCondition
{
public:
   MFEMInitialConditionFromFolder(FEModel* fem) : MFEMInitialCondition(fem) {}

   void Apply(MFEMSolver* solver) override;

private:
   std::string m_folder_path = "InitialConditions";

   DECLARE_FECORE_CLASS();
};
