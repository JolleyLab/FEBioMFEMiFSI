#pragma once
#include <FEBioMech/FEElasticMaterialPoint.h>

class MFEMImmersedElasticMaterialPoint : public FEElasticMaterialPoint
{
public:
    //! constructor
    MFEMImmersedElasticMaterialPoint(FEMaterialPointData* pt = 0);
    
public:
    // fluid data
    mat3d gradv; // spatial velocity gradient
    mat3ds ev; // strain rate
};
