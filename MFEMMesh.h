#pragma once
#include <mfem.hpp>
#include <map>

class FEMesh;
class FEFacetSet;

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

private:
	std::map<FEFacetSet*, int> m_surfMap;
};

MFEMMesh* CreateMFEMMesh3D(FEMesh& mesh);
