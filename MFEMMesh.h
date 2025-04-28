#pragma once
#include <mfem.hpp>
#include <map>

class FEMesh;
class FEFacetSet;
class FEDomain;
class FEModel;

class MFEMMesh : public mfem::Mesh 
{
public:
	MFEMMesh(int nodes, int elems, int bndrelems);

	void AssignSurfaceAttribute(FEFacetSet* set, int attr)
	{
		m_surfMap[set] = attr;
	}

	int GetSurfaceAttribute(FEFacetSet* set)
	{
		if (m_surfMap.find(set) == m_surfMap.end()) return -1;
		return m_surfMap[set];
	}

   void AssignDomainAttribute(FEDomain* dom, int attr)
   {
      m_domMap[dom] = attr;
   }

   int GetDomainAttribute(FEDomain* dom)
   {
      if (m_domMap.find(dom) == m_domMap.end()) return -1;
      return m_domMap[dom];
   }

private:
	std::map<FEFacetSet*, int> m_surfMap;
   std::map<FEDomain*, int> m_domMap;
};

MFEMMesh* CreateMFEMMesh3D(FEMesh& mesh);
