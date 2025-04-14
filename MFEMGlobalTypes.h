#pragma once

#include <FECore/FEGlobalVector.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/CompactUnSymmMatrix.h>

namespace mfem
{
   class SparseMatrix;
   class Vector;
   template <typename T> class Array;
}

class MFEMGlobalVector : public FEGlobalVector
{
public:
   // constructor
   MFEMGlobalVector(FEModel* fem, int size, std::vector<double>& R, std::vector<double>& Fr);
   // destructor
   ~MFEMGlobalVector() noexcept override;

   //! Assemble the element vector into this global vector
	void Assemble(std::vector<int>& en, std::vector<int>& elm, std::vector<double>& fe, bool bdom = false) override;

	//! Assemble into this global vector
	void Assemble(std::vector<int>& lm, std::vector<double>& fe) override;

	//! assemble a nodel value
	void Assemble(int node, int dof, double f) override;

   // Assemble into global mfem vector
   void Assemble(mfem::Array<int> &vdofs, std::vector<double>& fe);
   void Assemble(mfem::Array<int> &vdofs, mfem::Vector& fe);

   mfem::Vector &GetMfemVector() { return *mfem_R; }

public:
   mfem::Vector *mfem_R = nullptr;
};

class MFEMGlobalMatrix : public FEGlobalMatrix
{
public:
   MFEMGlobalMatrix(CRSSparseMatrix* pK, bool del = true);

   void CreateProfile(FEModel* pfem, int neq, bool reshape);
   void BuildMatrixProfile(FEModel* pfem, bool breset);

   mfem::SparseMatrix *GetMfemMatrix() { return mfem_K; }

   virtual ~MFEMGlobalMatrix() noexcept {};

private:
   mfem::SparseMatrix *mfem_K = nullptr;
};