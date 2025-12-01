#pragma once

#include "mfem.hpp"

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