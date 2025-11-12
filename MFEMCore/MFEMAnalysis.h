#pragma once
#include <FECore/FEAnalysis.h>

class MFEMAnalysis : public FEAnalysis
{
public:
	MFEMAnalysis(FEModel* fem);
	DECLARE_FECORE_CLASS();
};
