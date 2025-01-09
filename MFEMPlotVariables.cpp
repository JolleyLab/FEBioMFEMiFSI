#include "MFEMPlotVariables.h"
#include <FECore/FEModel.h>

bool FEPltMFemPressure::Save(FEMesh& m, FEDataStream& a)
{
	int dofP = GetFEModel()->GetDOFIndex("p");
	for (int i = 0; i < m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		a << node.get(dofP);
	}
	return true;
};

bool FEPltMFemVelocity::Save(FEMesh& m, FEDataStream& a)
{
	int dofVX = GetFEModel()->GetDOFIndex("vx");
	int dofVY = GetFEModel()->GetDOFIndex("vy");
	int dofVZ = GetFEModel()->GetDOFIndex("vz");
	for (int i = 0; i < m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		a << node.get(dofVX) << node.get(dofVY) << node.get(dofVZ);
	}
	return true;
};
