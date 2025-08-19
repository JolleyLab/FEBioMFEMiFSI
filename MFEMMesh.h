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

   void SetFEBioElementID(int mfem_id, int febio_id)
   {
      m_mfem_2_febio_elemMap[mfem_id] = febio_id;
   };

   int GetFEBioElementID(int mfem_id)
   {
      if (m_mfem_2_febio_elemMap.find(mfem_id) == m_mfem_2_febio_elemMap.end()) return -1;
      return m_mfem_2_febio_elemMap[mfem_id];
   };

   void SetNumSurfaces(int n)
   {
      m_mfem_surfelem_to_febio_surfelemMap.resize(n);
   };

   void SetFEBioSurfaceElementID(int surf_attr, int mfem_id, int febio_id)
   {
      m_mfem_surfelem_to_febio_surfelemMap[surf_attr-1][mfem_id] = febio_id;
   };

   int GetFEBioSurfaceElementID(int surf_attr, int mfem_id)
   {
      if (surf_attr-1 < 0 || surf_attr-1 >= (int)m_mfem_surfelem_to_febio_surfelemMap.size()) return -2;
      if (m_mfem_surfelem_to_febio_surfelemMap[surf_attr-1].find(mfem_id) == m_mfem_surfelem_to_febio_surfelemMap[surf_attr-1].end()) return -1;
      return m_mfem_surfelem_to_febio_surfelemMap[surf_attr-1][mfem_id];
   };

private:
	std::map<FEFacetSet*, int> m_surfMap;
   std::map<FEDomain*, int> m_domMap;
   std::map<int,int> m_mfem_2_febio_elemMap;
   std::vector<std::map<int,int>> m_mfem_surfelem_to_febio_surfelemMap;
};

MFEMMesh* CreateMFEMMesh3D(FEMesh& mesh);
