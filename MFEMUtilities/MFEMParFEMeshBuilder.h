#pragma once

#include <FECore/FEModel.h>

class MFEMMesh;

#include <mfem.hpp>

#ifdef MFEM_USE_MPI
class MFEMParFEMeshBuilder
{
public:
   MFEMParFEMeshBuilder(FEModel* fem, MFEMMesh* mesh_);
   ~MFEMParFEMeshBuilder();

   void SetMFEMSubMesh(mfem::SubMesh* submesh)
   {
      m_submesh = submesh;
   };

   void BuildLocalFEBioMesh(mfem::ParFiniteElementSpace* dfes);

   // Updates domain data to contain only the elements within a partition
   // specified by mfem ParMesh
   void ParCopyDomain(FEDomain &dom, FEDomain &par_dom, mfem::ParFiniteElementSpace &pfes);

   void ClearSerialSurfaces();

private:
   // FE model
   FEModel *m_fem = nullptr;

   // Global MFEM serial mesh
   MFEMMesh* m_serial_mesh = nullptr;

   // Global FEBio serial mesh
   FEMesh* m_serial_femesh = nullptr;

   // Serial submesh
   mfem::SubMesh* m_submesh = nullptr;

   // FEBio local surface mesh
   FESurface* m_febio_surf = nullptr;

   // Array of FESurfaces
   mfem::Array<FESurface*> serial_surf_list;

   // Array of maps from FEBio element id to local mfem element id
   std::vector<mfem::Array<int>> m_list_elem_id_maps;
};
#endif