#include "MFEMImmersedElasticMaterialPoint.h"

MFEMImmersedElasticMaterialPoint::MFEMImmersedElasticMaterialPoint(FEMaterialPointData* pt) : FEElasticMaterialPoint(pt)
{
   gradv.zero();
   ev.zero();
}
