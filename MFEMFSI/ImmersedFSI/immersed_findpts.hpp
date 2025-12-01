#pragma once

#include "mfem.hpp"

namespace immersed_solid
{
#ifdef MFEM_USE_GSLIB
// special findpts to do force spreading for immersed fsi applications, assumes every processor has part of fluid and solid mesh
class ImmersedFSI_FindPts : public mfem::FindPointsGSLIB
{
private:
   #ifdef MFEM_USE_MPI
   HYPRE_ParCSRMatrix pcsr_matrix;
   HYPRE_IJMatrix ij_matrix;
   mfem::HypreParMatrix *C = nullptr;
   mfem::HypreParMatrix *CT = nullptr;

   // Map from local dof of the search mesh to the master rank
   mfem::ParGridFunction my_rank_gf;
   // Map local tdof on solid mesh to ldof
   std::map<int,int> ltdof_to_ldof;
   // max rowsize for ij matrix
   int max_rowsize[1];
   #else
   mfem::SparseMatrix *C = nullptr;
   mfem::SparseMatrix *CT = nullptr;
   #endif

public:
   ImmersedFSI_FindPts();

#ifdef MFEM_USE_MPI
   ImmersedFSI_FindPts(MPI_Comm comm_);
#endif

   /** Spreads the fsi force from solid nodes to fluid nodes 
    * field_in = solid nodal forces
    * field_out = fluid nodal forces 
    * point_pos = solid nodes
    * 
    * the output field will need a group summation to work in parallel
   */
   void SpreadForces(const mfem::Vector &point_pos, const mfem::Vector &field_in, mfem::GridFunction &field_out, int point_pos_ordering = mfem::Ordering::byNODES);
   void Spread(const mfem::Vector &field_in, mfem::GridFunction &field_out);

   void TransferIndicator(const mfem::Vector &point_pos, mfem::GridFunction &field_out, int point_pos_ordering = mfem::Ordering::byNODES);
   void Transfer(mfem::GridFunction &field_out);

   void AssembleCouplingMatrix(const mfem::Vector &solid_nodes, mfem::FiniteElementSpace *fes_in, mfem::FiniteElementSpace *fes_out, int orderingtype_);
   #ifdef MFEM_USE_MPI
   void AssembleC(mfem::FiniteElementSpace *fes_in, mfem::FiniteElementSpace *fes_out, const mfem::Vector &solid_nodes, int orderingtype);
   #else
   void AssembleC(mfem::FiniteElementSpace *fes_in, mfem::FiniteElementSpace *fes_out);
   #endif

   mfem::Operator &GetCouplingMatrix()
   {
      return *C;
   };

   mfem::Operator &GetCouplingMatrixTranspose()
   {
      return *CT;
   };
   
   #ifdef MFEM_USE_MPI
   // Creates a map from a given dof of the search mesh to the master rank
   void SetupLocalDofToMasterRankMap(mfem::ParFiniteElementSpace *fes);

   // Determine the max row size for C matrix (i.e. max number of fluid dofs a solid node may influence)
   void DetermineMaxRowSize(mfem::ParFiniteElementSpace *fes);

   // Creates a map from a local tdof in the solid to its local dof number
   void SetupLocalTDofToLDofSolidMap(mfem::ParFiniteElementSpace *solid_fes);
   #endif

   // destructor
   virtual ~ImmersedFSI_FindPts();
};
#endif
} // namespace immersed_solid





