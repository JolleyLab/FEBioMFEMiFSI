#pragma once

#include "mfem.hpp"


#ifdef MFEM_USE_MPI
// extended ParAdd function to accept coefficients when adding two matrices
mfem::HypreParMatrix * ParAddCoeffs(double alpha, const mfem::HypreParMatrix *A, double beta, const mfem::HypreParMatrix *B);
void HypreParMatrixScaleRows(mfem::HypreParMatrix *B, const mfem::Vector &diag, const double &alpha);
void HypreParMatrixInvScaleRows(mfem::HypreParMatrix *B, const mfem::Vector &diag, const double &alpha);
void HypreParMatrixAddToDiag(mfem::HypreParMatrix *B, const mfem::Vector &diag);
#endif
