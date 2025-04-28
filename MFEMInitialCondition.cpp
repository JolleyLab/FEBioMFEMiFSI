#include "MFEMInitialCondition.h"
#include "MFEMSolver.h"
#ifdef WIN32
#include <mfem/fem/gridfunc.hpp>
#else
#include "mfem.hpp"
#endif
#include "MFEMTools.h"

BEGIN_FECORE_CLASS(MFEMInitialPressure, MFEMInitialCondition)
	ADD_PARAMETER(m_pressure, "pressure");
END_FECORE_CLASS();

void MFEMInitialPressure::Apply(MFEMSolver* solver)
{
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction* p_ic = solver->GetCurrentPressure();
   #else
	mfem::GridFunction* p_ic = solver->GetCurrentPressure();
   #endif
	FEScalarCoefficient p_coeff(GetFEModel(), m_pressure);
	p_coeff.SetTime(0.0);
	p_ic->ProjectCoefficient(p_coeff);
}

void MFEMInitialPressure::CopyFrom(FEInitialCondition& ic)
{
   MFEMInitialCondition::CopyFrom(ic);

   MFEMInitialPressure& p = dynamic_cast<MFEMInitialPressure&>(ic);
   CopyParameterListState(p.GetParameterList());

   (*this).m_pressure = p.m_pressure;
}

double MFEMInitialPressure::GetPressure() const
{ 
   if (m_pressure.isConst()) return m_pressure.constValue();
   else
   {
      MFEM_ABORT("MFEMInitialPressure::GetPressure: non-constant pressure");
      return 0.0;
   }
}

BEGIN_FECORE_CLASS(MFEMInitialVelocity, MFEMInitialCondition)
	ADD_PARAMETER(m_vx, "vx");
	ADD_PARAMETER(m_vy, "vy");
	ADD_PARAMETER(m_vz, "vz");
END_FECORE_CLASS();

void MFEMInitialVelocity::Apply(MFEMSolver* solver)
{
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction* u_ic = solver->GetCurrentVelocity();
   #else
	mfem::GridFunction* u_ic = solver->GetCurrentVelocity();
   #endif
	FEVectorCoefficient u_coeff(GetFEModel(), m_vx, m_vy, m_vz);
	u_coeff.SetTime(0.0);
	u_ic->ProjectCoefficient(u_coeff);
}

void MFEMInitialVelocity::CopyFrom(FEInitialCondition& ic)
{
   MFEMInitialCondition::CopyFrom(ic);

   MFEMInitialVelocity& v = dynamic_cast<MFEMInitialVelocity&>(ic);
   CopyParameterListState(v.GetParameterList());
   
   (*this).m_vx = v.m_vx;
   (*this).m_vy = v.m_vy;
   (*this).m_vz = v.m_vz;
}
