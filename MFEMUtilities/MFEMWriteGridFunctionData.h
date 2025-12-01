#pragma once

#include <FECore/FEElement.h>
#include <FECore/FEMaterialPoint.h>
#include <functional>

#include "mfem.hpp"

template <class T> void computeNodalProjectedValues(FEElement& e, mfem::GridFunction& gf, mfem::Array<int>& zones_per_vdof, mfem::Array<int>& vdofs, std::function<T(const FEMaterialPoint&)> var)
{
	T gi[FEElement::MAX_INTPOINTS];
	T gn[FEElement::MAX_NODES];

	T si[FEElement::MAX_INTPOINTS];
	T sn[FEElement::MAX_NODES];

   int ne = e.Nodes();
   int ni = e.GaussPoints();
   int ndofs = vdofs.Size();

   // get the integration point values
   for (int k = 0; k<ni; ++k)
   {
      FEMaterialPoint& mp = *e.GetMaterialPoint(k);
      T s = var(mp);
      si[k] = s;
   }

   // project to nodes
   e.project_to_nodes(si, sn);

   // add data to gf and count zones
   for (int j=0; j<ndofs; ++j)
   {
      gf(vdofs[j]) += sn[j]; 
      zones_per_vdof[vdofs[j]]++;
   }     
}