#pragma once

#include "mfem.hpp"

namespace ins_solver
{
class BlockOrthoSolver : public mfem::Solver
{
private:
   #ifdef MFEM_USE_MPI
   MPI_Comm mycomm;
   mutable HYPRE_BigInt global_size;
   #else
   mutable int global_size;
   #endif
   const bool parallel;
   mfem::Array<int> bOffsets, vblock, pblock;

public:
   BlockOrthoSolver(mfem::Array<int> &bOffsets_);
   #ifdef MFEM_USE_MPI
   BlockOrthoSolver(mfem::Array<int> &bOffsets_, MPI_Comm mycomm_);
   #endif

   void SetSolver(mfem::Solver &s);

   mfem::Solver *GetSolver() { return solver; };

   virtual void SetOperator(const mfem::Operator &op);

   void Mult(const mfem::Vector &b, mfem::Vector &x) const;

   int GetNumIterations()
   {
      return dynamic_cast<mfem::IterativeSolver*>(solver)->GetNumIterations();
   };
   double GetFinalNorm()
   {
      return dynamic_cast<mfem::IterativeSolver*>(solver)->GetFinalNorm();
   };
   double GetFinalRelNorm()
   {
      return dynamic_cast<mfem::IterativeSolver*>(solver)->GetFinalRelNorm();
   }; 

private:
	mfem::Solver *solver = nullptr;

   mutable mfem::Vector b_ortho, p_ortho, p_vec, temp;

   void Orthogonalize(const mfem::Vector &v, mfem::Vector &v_ortho) const;
};

}
