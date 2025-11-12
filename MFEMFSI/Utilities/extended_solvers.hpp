#pragma once

#include "mfem.hpp"

#ifdef MFEM_USE_MPI
#include "../Elasticity/SpatialOperators/residual_op.hpp"
#include "../Elasticity/bc_containers.hpp"
#endif

// class to extend MFEM newton solver with additional options
class ExtendedNewtonSolver : public mfem::NewtonSolver
{
protected:
   mutable bool apply_min_iter_after1 = false;
   int min_iter = 0; // min number of iterations to do (0 by default)
   bool gen_alpha = false;
   mfem::Vector *u_np1 = nullptr;
   double gamma;
   double *dt = nullptr;
   mutable mfem::BlockVector kblock;
   mfem::Array<int> block_trueOffsets;

public:
   ExtendedNewtonSolver() : NewtonSolver() {};

   #ifdef MFEM_USE_MPI
   ExtendedNewtonSolver(MPI_Comm comm_) : NewtonSolver(comm_) {};
   #endif

   // extended mult method
   virtual void Mult(const mfem::Vector &b, mfem::Vector &x) const;

   void SetMinIter(int min_it)
   {
      min_iter = min_it;
   };

   void ApplyMinIterAfterFirstStep()
   {
      apply_min_iter_after1 = true;
   };

   void SetGenAlphaUpdate(mfem::Vector &un_next_, double gamma_, double &curr_dt, mfem::Array<int> &block_trueOffsets_)
   {
      gen_alpha = true;
      u_np1 = &un_next_;
      gamma = gamma_;
      dt = &curr_dt;
      block_trueOffsets = block_trueOffsets_;
   }

   void UpdateGenAlphaVelocity(const mfem::Vector &k) const
   {
      kblock.Update(k.GetData(),block_trueOffsets); // delta [udot p]
      const mfem::Vector &delta_udot = kblock.GetBlock(0);

      // u_n+1,(i) = u_n+1, (i-1) + gamma*dt*DELTA_udot_n+1,(i)
      u_np1->Add(gamma*(*dt),delta_udot);
   }

   virtual ~ExtendedNewtonSolver() {};
};

class ExtendedFGMRESSolver : public mfem::FGMRESSolver
{
protected:
   int min_iter = 1; // min number of iterations to do (1 by default)

public:
   ExtendedFGMRESSolver() : mfem::FGMRESSolver() {};
   
   #ifdef MFEM_USE_MPI
   ExtendedFGMRESSolver(MPI_Comm comm_) : mfem::FGMRESSolver(comm_) { };
   #endif

   // extended mult method
   virtual void Mult(const mfem::Vector &b, mfem::Vector &x) const;

   void SetMinIter(int min_it)
   {
      min_iter = min_it;
   }
};

class ExtendedGMRESSolver : public mfem::GMRESSolver
{
protected:
   int min_iter = 0; // min number of iterations to do (0 by default)

public:
   ExtendedGMRESSolver() : mfem::GMRESSolver() {};

   #ifdef MFEM_USE_MPI
   ExtendedGMRESSolver(MPI_Comm comm_) : GMRESSolver(comm_) { };
   #endif

   // extended mult method
   virtual void Mult(const mfem::Vector &b, mfem::Vector &x) const;

   void SetMinIter(int min_it)
   {
      min_iter = min_it;
   }
};

/// Modified GMRES method
// allows for preconditioner to use a different operator for the preconditioner
// and the iterative solver. This is useful when it is more efficient to utilize a matrix-free
// version of the original operator and a sparse approx of the operator for the preconditioner
namespace modified_gmres_serial
{
   class ModifiedGMRESSolver : public mfem::GMRESSolver
   {
   protected:
      int m; // see SetKDim()

   public:
      ModifiedGMRESSolver() : mfem::GMRESSolver() { m = 50; };

   #ifdef MFEM_USE_MPI
      ModifiedGMRESSolver(MPI_Comm comm_) : mfem::GMRESSolver(comm_) { m = 50; };
   #endif

      /// Set the number of iteration to perform between restarts, default is 50.
      void SetKDim(int dim) { m = dim; }

      virtual void Mult(const mfem::Vector &b, mfem::Vector &x) const;

      void SetOperatorGMRES(const mfem::Operator &op_gmres_)
      {
         oper = &op_gmres_;
         height = op_gmres_.Height();
         width = op_gmres_.Width();
      };
      void SetOperatorPrec(const mfem::Operator &op)
      {
         prec->SetOperator(op);
      };
   };
}