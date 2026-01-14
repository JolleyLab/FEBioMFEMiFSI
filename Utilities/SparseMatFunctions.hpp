#pragma once

#include "mfem.hpp"

void SparseMatInvScaleRows(mfem::SparseMatrix & B, const mfem::Vector & sl, double alpha);
void SparseMatScaleRows(mfem::SparseMatrix & B, const mfem::Vector & sl, double alpha);
void SparseMatAddToDiag(mfem::SparseMatrix & B, const mfem::Vector & d);