#pragma once
#include <FECore/FEMaterial.h>

class MFEMMaterial : public FEMaterial
{
public:
	MFEMMaterial(FEModel* fem);

	DECLARE_FECORE_CLASS();
};
