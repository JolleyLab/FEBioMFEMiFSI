#pragma once
#include <FECore/FEMaterial.h>

class MFEMMaterial : public FEMaterial
{
public:
	MFEMMaterial(FEModel* fem);

   // return material properties
   virtual double Density() = 0;
};