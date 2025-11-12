#pragma once

#include "mfem.hpp"

#ifdef MFEM_USE_ZLIB
class SimvascularMeshReader : public mfem::Mesh
{
public:
   SimvascularMeshReader(const std::string &filename, int generate_edges = 0, int refine = 1, bool fix_orientation = true);

   void Load(std::istream &input, int generate_edges = 0, int refine = 1, bool fix_orientation = true) override;

   void SVLoader(std::istream &input, int generate_edges = 0);

   void SV_ReadXML_VTKMesh(std::istream &input, int &curved, int &read_gf,
      bool &finalize_topo, const std::string &xml_prefix);
};
#endif