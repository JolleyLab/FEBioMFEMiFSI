#include "HypreFunctions.hpp"
#include "mfem.hpp"
#include <general/forall.hpp>

using namespace std;
using namespace mfem;

#ifdef MFEM_USE_MPI
HypreParMatrix * ParAddCoeffs(double alpha, const HypreParMatrix *A, double beta, const HypreParMatrix *B)
{
   hypre_ParCSRMatrix *C;
#if MFEM_HYPRE_VERSION <= 22000
   hypre_ParcsrAdd(alpha, *A, beta, *B, &C);
#else
   hypre_ParCSRMatrixAdd(alpha, *A, beta, *B, &C);
#endif
   if (!hypre_ParCSRMatrixCommPkg(C)) { hypre_MatvecCommPkgCreate(C); }

   return new HypreParMatrix(C);
}

void HypreParMatrixInvScaleRows(HypreParMatrix *B, const Vector &diag, const double &alpha)
{
   hypre_ParCSRMatrix *A = static_cast<hypre_ParCSRMatrix*>(*B);
   if (hypre_CSRMatrixNumRows(A->diag) != hypre_CSRMatrixNumRows(A->offd))
   {
      mfem_error("Row does not match");
   }

   if (hypre_CSRMatrixNumRows(A->diag) != diag.Size())
   {
      mfem_error("Note the Vector diag is not of compatible dimensions with A\n");
   }

   B->HostReadWrite();
   diag.HostRead();

   int size = B->Height();
   real_t     *Adiag_data   = hypre_CSRMatrixData(A->diag);
   HYPRE_Int  *Adiag_i      = hypre_CSRMatrixI(A->diag);


   real_t     *Aoffd_data   = hypre_CSRMatrixData(A->offd);
   HYPRE_Int  *Aoffd_i      = hypre_CSRMatrixI(A->offd);
   real_t val;
   HYPRE_Int jj;
   for (int i(0); i < size; ++i)
   {
#ifdef MFEM_DEBUG
      if (0.0 == diag(i))
      {
         mfem_error("HypreParMatrix::InvDiagScale : Division by 0");
      }
#endif
      val = 1./diag(i)/alpha;
      for (jj = Adiag_i[i]; jj < Adiag_i[i+1]; ++jj)
      {
         Adiag_data[jj] *= val;
      }
      for (jj = Aoffd_i[i]; jj < Aoffd_i[i+1]; ++jj)
      {
         Aoffd_data[jj] *= val;
      }
   }

   B->HypreRead();
}

void HypreParMatrixScaleRows(HypreParMatrix *B, const Vector &diag, const double &alpha)
{
   hypre_ParCSRMatrix *A = static_cast<hypre_ParCSRMatrix*>(*B);
   if (hypre_CSRMatrixNumRows(A->diag) != hypre_CSRMatrixNumRows(A->offd))
   {
      mfem_error("Row does not match");
   }

   if (hypre_CSRMatrixNumRows(A->diag) != diag.Size())
   {
      mfem_error("Note the Vector diag is not of compatible dimensions with A\n");
   }

   B->HostReadWrite();
   diag.HostRead();

   int size = B->Height();
   real_t     *Adiag_data   = hypre_CSRMatrixData(A->diag);
   HYPRE_Int  *Adiag_i      = hypre_CSRMatrixI(A->diag);

   real_t     *Aoffd_data   = hypre_CSRMatrixData(A->offd);
   HYPRE_Int  *Aoffd_i      = hypre_CSRMatrixI(A->offd);
   real_t val;
   HYPRE_Int jj;
   for (int i(0); i < size; ++i)
   {
      val = diag[i] * alpha;
      for (jj = Adiag_i[i]; jj < Adiag_i[i+1]; ++jj)
      {
         Adiag_data[jj] *= val;
      }
      for (jj = Aoffd_i[i]; jj < Aoffd_i[i+1]; ++jj)
      {
         Aoffd_data[jj] *= val;
      }
   }

   B->HypreRead();
}

void HypreParMatrixAddToDiag(HypreParMatrix *B, const Vector &diag)
{
   hypre_ParCSRMatrix *A = static_cast<hypre_ParCSRMatrix*>(*B);
   
   const int size = diag.Size();
   auto hypre_ml = GetHypreMemoryLocation();
   // Avoid using GetHypreMemoryClass() since it may be MemoryClass::MANAGED and
   // that may not play well with the memory types used by 'diag'.
   MemoryClass hypre_mc = (hypre_ml == HYPRE_MEMORY_HOST) ?
                          MemoryClass::HOST : MemoryClass::DEVICE;
   const real_t *diag_hd = diag.GetMemory().Read(hypre_mc, size);
#if MFEM_HYPRE_VERSION >= 21800
   MFEM_VERIFY(A->diag->memory_location == hypre_ml,
               "unexpected HypreParMatrix memory location!");
#endif
   const HYPRE_Int *A_diag_i = A->diag->i;
   real_t *A_diag_d = A->diag->data;
#ifdef MFEM_DEBUG
   const HYPRE_Int *A_diag_j = A->diag->j;
#endif
   mfem::hypre_forall(size, [=] MFEM_HOST_DEVICE (int i)
   {
      A_diag_d[A_diag_i[i]] += diag_hd[i];
      MFEM_ASSERT_KERNEL(
         A_diag_j[A_diag_i[i]] == i,
         "The first entry in each row must be the diagonal one!");
   });
}
#endif