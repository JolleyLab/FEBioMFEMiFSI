#pragma once

// residual operator
#include "mfem.hpp"
#include "CVBC_AinvPC.hpp"
#include "../../../Utilities/hypre_func_extensions.hpp"
#include "../../../Utilities/extended_solvers.hpp"

// Base Class for Implementing 2x2 Block Preconditioners for matrices of the form
//        Matrix System = [ A   G ]
//                        [ D   C ]
// Here we assume an LDU factorization
class BlockPreconditioner : public mfem::Solver
{
protected:
   // Finite element spaces for setting up preconditioner blocks
   mfem::Array<mfem::FiniteElementSpace *> spaces;
   
   // temp var to apply ops
   mutable mfem::BlockVector xblock;
   mutable mfem::BlockVector yblock;
   mutable mfem::Vector tmp1,tmp2,tmp3;
   #ifdef MFEM_USE_MPI
   mfem::HypreParMatrix *tempG = nullptr;
   mfem::HypreParMatrix *tempS = nullptr; // owned
   mfem::HypreParMatrix *tempA = nullptr;
   mfem::HypreParMatrix *tempD = nullptr;
   mfem::HypreParMatrix *tempC = nullptr;
   mfem::HypreParMatrix *S = nullptr; // schur complement operator // owned
   #else
   mfem::SparseMatrix *tempG = nullptr;
   mfem::SparseMatrix *tempS = nullptr; // owned
   mfem::SparseMatrix *tempA = nullptr;
   mfem::SparseMatrix *tempD = nullptr;
   mfem::SparseMatrix *tempC = nullptr;
   mfem::SparseMatrix *S = nullptr; // schur complement operator // owned
   #endif
   mfem::Array<int> offsets;
   double Gcoeff = 1.0;
   double Dcoeff = 1.0;

   // vector diagonal of A operator
   mfem::Vector Ad, Asparse;

   // Jacobian for block access
   mfem::BlockOperator *jacobian = nullptr;

   // Preconditioners and Solvers
   mfem::Solver *AInv = nullptr; // owned
   mfem::Solver *AInvPC = nullptr; // owned
   mfem::Solver *SInv = nullptr; // owned
   mfem::Solver *SInvPC = nullptr; // owned

   // options
   bool equal_order = false;
   #ifdef MFEM_USE_MPI
   bool jacobi_preconditioning = false; // default is AMG
   #else
   bool jacobi_preconditioning = true; // default is Jacobi in serial
   #endif
   bool direct_solvers = false; // use direct solvers for A and S blocks

   double a_atol = 1E-10;
   double a_rtol = 1E-6;
   int a_maxIt = 200;

   double s_atol = 1E-10;
   double s_rtol = 1E-6;
   int s_maxIt = 200;

   // CV bcs 
   bool cv_bcs = false;
   CVBC_AinvPC *cvbc_PC = nullptr; // not owned
   bool svsolver = true; // use AinvPC from Simvascular solver (https://doi.org/10.1007/s00466-013-0868-1)
   #ifdef MFEM_USE_MPI
   mutable mfem::HypreParMatrix *A_cvbc = nullptr; // owned 
   #else
   mutable mfem::SparseMatrix *A_cvbc = nullptr; // owned
   #endif
public:
   BlockPreconditioner(mfem::Array<int> &offsets_, mfem::Array<mfem::FiniteElementSpace *> &fes, bool equal_order_);

   // Setup functions
   virtual void SetOperator(const mfem::Operator &op);
   virtual void SetASolver(mfem::BlockOperator *jacobian);
   virtual void SetSchurSolver(mfem::BlockOperator *jacobian);
   virtual void FormSchurOperator(mfem::BlockOperator *jacobian) = 0;

   // note this is only used when there are no CVBCs
   void SetAinvTol(double a_atol_, double a_rtol_, int a_maxIt_ = 200)
   {
      a_atol = a_atol_;
      a_rtol = a_rtol_;
      a_maxIt = a_maxIt_;
   };

   void SetSinvTol(double s_atol_, double s_rtol_, int s_maxIt_ = 200)
   {
      s_atol = s_atol_;
      s_rtol = s_rtol_;
      s_maxIt = s_maxIt_;
   };

   void SetCV_Ainv_PC(CVBC_AinvPC *cvbc_PC_)
   {
      cvbc_PC = cvbc_PC_;
      cv_bcs = true;
   };

   /// Operator application
   virtual void Mult (const mfem::Vector & x, mfem::Vector & y) const = 0;

   /// Action of the transpose operator
   virtual void MultTranspose (const mfem::Vector & x, mfem::Vector & y) const;

   void UseJacobiPreconditioning()
   {
      jacobi_preconditioning = true;
   }

   void UseDirectSolvers()
   {
      #ifndef MFEM_USE_SUITESPARSE
      MFEM_ABORT("Error: Direct solvers require mfem with SuiteSparse!");
      #endif
      direct_solvers = true;
   }

   virtual ~BlockPreconditioner();
};