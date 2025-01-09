#pragma once
#include <FECore/FEPlotData.h>

class FEPltMFemPressure : public FEPlotNodeData
{
public:
	FEPltMFemPressure(FEModel* fem) : FEPlotNodeData(fem, PLT_FLOAT, FMT_NODE) {}
	bool Save(FEMesh& m, FEDataStream& a) override;
};

class FEPltMFemVelocity : public FEPlotNodeData
{
public:
	FEPltMFemVelocity(FEModel* fem) : FEPlotNodeData(fem, PLT_VEC3F, FMT_NODE) {}
	bool Save(FEMesh& m, FEDataStream& a) override;
};

