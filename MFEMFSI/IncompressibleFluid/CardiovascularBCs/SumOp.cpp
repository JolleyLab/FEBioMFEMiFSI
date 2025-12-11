#include "SumOp.hpp"
#include <general/forall.hpp>

using namespace mfem;

void SumOp::Mult(const Vector &k, Vector &y) const
{
   A->Mult(k,y);
   B->AddMultGrad(k,y,cv_time_factor);
   if (C)
   {
      C->AddMult(k,y);
   }

   // copy essential dofs in k to y
   // velocity
   const int N = vel_ess_tdofs->Size();
   const auto idx = vel_ess_tdofs->Read();
   auto y_RW = y.ReadWrite();
   auto k_R = k.Read();
   MFEM_FORALL(i, N, y_RW[idx[i]] = k_R[idx[i]];);
}