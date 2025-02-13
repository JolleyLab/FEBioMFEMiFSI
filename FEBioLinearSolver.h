#pragma once

#include <mfem/linalg/operator.hpp>
#include <mfem/linalg/sparsemat.hpp>
#include <FECore/LinearSolver.h>
#include <FECore/CompactUnSymmMatrix.h>

// an adaptor for using an mfem::SparseMatrix as an (FEBio) SparseMatrix
class MFEMSparseMatrix : public CRSSparseMatrix
{
public:
	MFEMSparseMatrix(mfem::SparseMatrix* A);
	~MFEMSparseMatrix();

private:
	mfem::SparseMatrix* m_mat;
};

// Adapter that allows the fluid solver to use the FEBio linear solvers
class FEBioLinearSolver : public mfem::Solver
{
public:
	FEBioLinearSolver(FEModel* fem);

	void SetOperator(const mfem::Operator& op) override;

	void Mult(const mfem::Vector& x, mfem::Vector& y) const override;

private:
	MFEMSparseMatrix* mat;
	LinearSolver* linSolver;
};
