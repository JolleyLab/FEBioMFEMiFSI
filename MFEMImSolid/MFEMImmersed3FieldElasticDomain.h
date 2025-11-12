#pragma once
#include "MFEMImmersedElasticDomain.h"

class MFEMImmersed3FieldElasticDomain : public MFEMImmersedElasticDomain
{
protected:
	struct ELEM_DATA
	{
		double	eJ;		// average element jacobian at intermediate time
		double	ep;		// average pressure
		double	Lk;		// Lagrangian multiplier
      double  eJt;    // average element jacobian at current time
      double  eJp;    // average element jacobian at previous time
	};

public:
   // constructor
   MFEMImmersed3FieldElasticDomain(FEModel* pfem);

   // destructor
   virtual ~MFEMImmersed3FieldElasticDomain() noexcept;

   //! initialize class
	bool Init() override;

   //! Reset data
	void Reset() override;

   //! initialize elements
   void PreSolveUpdate(const FETimeInfo& timeInfo) override;

public: // overriden from MFEMImmersedElasticDomain

   // solid stiffness
   void ElementSolidStiffness(int iel, FESolidElement& el, matrix& ke, double scale) override;

   //! geometrical stiffness (i.e. initial stress)
	void ElementGeometricalStiffness(int iel, FESolidElement& el, matrix& ke, double scale) override;

	//! material stiffness component
	void ElementMaterialStiffness(int iel, FESolidElement& el, matrix& ke, double scale) override;

   //! update the stress of an element
	void UpdateElementStress(int iel, const FETimeInfo& tp) override;

protected:

   //! Dilatational stiffness component for nearly-incompressible materials
	void ElementDilatationalStiffness(FEModel& fem, int iel, matrix& ke, double scale);


protected:
   vector<ELEM_DATA>	m_Data;
};


