#include "FEBioLinearSolver.h"
#include <FECore/FECoreKernel.h>

MFEMSparseMatrix::MFEMSparseMatrix(mfem::SparseMatrix* A) : m_mat(A) 
{
	// I think most FEBio solvers will expect indices to be sorted
	A->SortColumnIndices();

	m_offset = 0;
	m_nrow = A->Width();
	m_ncol = A->Height();
	m_nsize = A->NumNonZeroElems();

	const int* Ap = A->HostReadI();
	const int* Ai = A->HostReadJ();
	const mfem::real_t* Ax = A->HostReadData();

	m_ppointers = const_cast<int*>(Ap);
	m_pindices = const_cast<int*>(Ai);
	m_pd = const_cast<double*>(Ax);
}

MFEMSparseMatrix::~MFEMSparseMatrix()
{
	m_pindices = nullptr;
	m_ppointers = nullptr;
	m_pd = nullptr;
}


FEBioLinearSolver::FEBioLinearSolver(FEModel* fem) 
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	linSolver = fecore.CreateDefaultLinearSolver(fem);
	mat = nullptr;
}

void FEBioLinearSolver::SetOperator(const mfem::Operator& op)
{
	mfem::SparseMatrix* A = const_cast<mfem::SparseMatrix*>(dynamic_cast<const mfem::SparseMatrix*>(&op));
	assert(A);

	linSolver->Destroy();
	if (mat) delete mat;

	mat = new MFEMSparseMatrix(A);
	bool b = linSolver->SetSparseMatrix(mat); assert(b);

	b = linSolver->PreProcess(); assert(b);
	b = linSolver->Factor(); assert(b);
}

void FEBioLinearSolver::Mult(const mfem::Vector& b, mfem::Vector& x) const
{
	b.HostRead();
	x.HostReadWrite();
	int n = b.Size();
	assert(b.Size() == n);
	const mfem::real_t* pb = b.HostRead();
	std::vector<double> X(n), B(n);
	for (int i = 0; i < n; ++i) B[i] = pb[i];
	linSolver->BackSolve(X, B);
	mfem::real_t* px = x.HostWrite();
	for (int i = 0; i < n; ++i) px[i] = X[i];
}
