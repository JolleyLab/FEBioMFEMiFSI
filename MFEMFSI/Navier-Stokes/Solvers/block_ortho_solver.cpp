// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "block_ortho_solver.hpp"
#include <general/forall.hpp>

using namespace mfem;
using namespace ins_solver;

BlockOrthoSolver::BlockOrthoSolver(Array<int> &bOffsets_)
   : Solver(0, false), bOffsets(bOffsets_), global_size(-1), parallel(false)
{ 
   // form index array corresponding to velocity block
   int start_ind = bOffsets[0];
   int vblock_size = bOffsets[1] - bOffsets[0];
   vblock.SetSize(vblock_size);
   for (int i=0; i < vblock_size; i++)
   {
      vblock[i] = start_ind + i;
   }
   
   // form index array corresponding to pressure block
   start_ind = bOffsets[1];
   int pblock_size = bOffsets[2] - bOffsets[1];
   pblock.SetSize(pblock_size);
   for (int i=0; i < pblock_size; i++)
   {
      pblock[i] = start_ind + i;
   }

   p_vec.SetSize(pblock.Size());
   p_vec = 0.0;
   p_ortho = p_vec;

   // size temp vector
   temp.SetSize(vblock_size + pblock_size);
   temp = 0.0;
}

#ifdef MFEM_USE_MPI
BlockOrthoSolver::BlockOrthoSolver(Array<int> &bOffsets_, MPI_Comm mycomm_)
   : Solver(0, false), bOffsets(bOffsets_), mycomm(mycomm_), global_size(-1), parallel(true) 
{
   // form index array corresponding to velocity block
   int start_ind = bOffsets[0];
   int vblock_size = bOffsets[1] - bOffsets[0];
   vblock.SetSize(vblock_size);
   for (int i=0; i < vblock_size; i++)
   {
      vblock[i] = start_ind + i;
   }
   
   // form index array corresponding to pressure block
   start_ind = bOffsets[1];
   int pblock_size = bOffsets[2] - bOffsets[1];
   pblock.SetSize(pblock_size);
   for (int i=0; i < pblock_size; i++)
   {
      pblock[i] = start_ind + i;
   }

   p_vec.SetSize(pblock.Size());
   p_vec = 0.0;
   p_ortho = p_vec;

   // size temp vector
   temp.SetSize(vblock_size + pblock_size);
   temp = 0.0;
}
#endif

void BlockOrthoSolver::SetSolver(Solver &s)
{
   solver = &s;
   height = s.Height();
   width = s.Width();
   MFEM_VERIFY(height == width, "Solver must be a square Operator!");
   global_size = -1; // lazy evaluated
}

void BlockOrthoSolver::SetOperator(const Operator &op)
{
   MFEM_VERIFY(solver, "Solver hasn't been set, call SetSolver() first.");
   solver->SetOperator(op);
   height = solver->Height();
   width = solver->Width();
   MFEM_VERIFY(height == width, "Solver must be a square Operator!");
   global_size = -1; // lazy evaluated
}

void BlockOrthoSolver::Mult(const Vector &b, Vector &x) const
{
   MFEM_VERIFY(solver, "Solver hasn't been set, call SetSolver() first.");
   MFEM_VERIFY(height == solver->Height(),
               "solver was modified externally! call SetSolver() again!");
   MFEM_VERIFY(height == b.Size(), "incompatible input Vector size!");
   MFEM_VERIFY(height == x.Size(), "incompatible output Vector size!");

   // Orthogonalize input pressure block
   b.GetSubVector(pblock,p_vec);
   Orthogonalize(p_vec, p_ortho);
   // temp = b;
   // temp.SetSubVector(pblock,p_ortho); // update pressure block

   // RB 1-15-24
   // copy over velocity block
   const int N = vblock.Size();
   const auto idx = vblock.Read();
   auto y_RW = temp.ReadWrite();
   auto k_R = b.Read();
   MFEM_FORALL(i, N, y_RW[idx[i]] = k_R[idx[i]];);
   
   // update pressure block
   temp.SetSubVector(pblock,p_ortho);    
   
   // Propagate iterative_mode to the solver:
   solver->iterative_mode = iterative_mode;

   // Apply the Solver
   solver->Mult(temp, x);

   // Orthogonalize output
   x.GetSubVector(pblock,p_vec);
   Orthogonalize(p_vec,p_ortho);
   x.SetSubVector(pblock,p_ortho);
}

void BlockOrthoSolver::Orthogonalize(const Vector &v, Vector &v_ortho) const
{
   if (global_size == -1)
   {
      global_size = height;
      #ifdef MFEM_USE_MPI
      if (parallel)
      {
         MPI_Allreduce(MPI_IN_PLACE, &global_size, 1, HYPRE_MPI_BIG_INT,
                     MPI_SUM, mycomm);
      }
      #endif
   }

   // TODO: GPU/device implementation

   double global_sum = v.Sum();

   #ifdef MFEM_USE_MPI
   if (parallel)
   {
      MPI_Allreduce(MPI_IN_PLACE, &global_sum, 1, MPI_DOUBLE, MPI_SUM, mycomm);
   }
   #endif

   double ratio = global_sum / static_cast<double>(global_size);
   v_ortho.SetSize(v.Size());
   v.HostRead();
   v_ortho.HostWrite();
   for (int i = 0; i < v_ortho.Size(); ++i)
   {
      v_ortho(i) = v(i) - ratio;
   }
}