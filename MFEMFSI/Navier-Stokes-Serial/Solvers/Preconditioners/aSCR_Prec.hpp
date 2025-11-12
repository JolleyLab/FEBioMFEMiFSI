#pragma once

// residual operator
#include "mfem.hpp"
#include "BlockPreconditioner.hpp"

using namespace std;

#ifdef MFEM_USE_MPI
class aSCR_Prec_Serial : public BlockPreconditioner
{
private:
   // temp var to apply ops
   mutable mfem::Vector bhatp,yhatu;

   // pressue underrelaxation parameter
   double alpha = 0.9;
   // double alpha = 1.0;

public:
   aSCR_Prec_Serial(mfem::Array<int> &offsets_,mfem::Array<mfem::FiniteElementSpace *> &fes, bool equal_order_);

   void FormSchurOperator(mfem::BlockOperator *jacobian) override;

   void Mult(const mfem::Vector & x, mfem::Vector & y) const override;

   virtual ~aSCR_Prec_Serial();
};
#endif