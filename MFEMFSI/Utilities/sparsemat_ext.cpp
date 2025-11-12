#include "sparsemat_ext.hpp"
#include <general/forall.hpp>

using namespace mfem;

void SparseMatInvScaleRows(SparseMatrix & B, const Vector & sl, double alpha)
{
   double scale;
   if (B.Finalized())
   {
      int j, end;
      auto I = B.GetI();
      auto A = B.GetData();

      for (int i=0; i < B.Height(); ++i)
      {
         end = I[i+1];
         scale = 1.0/sl(i) * alpha;
         for (j = I[i]; j < end; j++)
         {
            A[j] *= scale;
         }
      }
   }
   else
   {
      MFEM_ABORT("SparseMatrix must be finalized");
   }
};

void SparseMatScaleRows(SparseMatrix & B, const Vector & sl, double alpha)
{
   double scale;
   if (B.Finalized())
   {
      int j, end;
      auto I = B.GetI();
      auto A = B.GetData();

      for (int i=0; i < B.Height(); ++i)
      {
         end = I[i+1];
         scale = sl(i) * alpha;
         for (j = I[i]; j < end; j++)
         {
            A[j] *= scale;
         }
      }
   }
   else
   {
      MFEM_ABORT("SparseMatrix must be finalized");
   }
};

void SparseMatAddToDiag(SparseMatrix & B, const Vector & d)
{
   int height = B.Height();
   int width = B.Width();
   MFEM_VERIFY(height == width, "Matrix must be square, not height = "
               << height << ", width = " << width);
   MFEM_VERIFY(B.Finalized(), "Matrix must be finalized.");

   const auto II = B.ReadI();
   const auto JJ = B.ReadJ();
   auto AA = B.ReadWriteData();
   const auto dd = d.Read();

   mfem::forall(height, [=] MFEM_HOST_DEVICE (int i)
   {
      const int begin = II[i];
      const int end = II[i+1];
      int j;
      for (j = begin; j < end; j++)
      {
         if (JJ[j] == i)
         {
            AA[j] += dd[i];
            break;
         }
      }
   });
};