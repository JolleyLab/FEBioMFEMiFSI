#pragma once

// residual operator
#include "mfem.hpp"
#include "../../../Navier-Stokes/ifs_bc_containers.hpp"
#include <general/forall.hpp>

using namespace std;

// Implements the special preconditioner from simvascular for CV bcs from
// https://link.springer.com/article/10.1007/s00466-013-0868-1#page-1
// Assumes Atotal = A + Abc where A is from 3D domain and Abc are the CV bc rank one matrices
//
// Future work could look at better than Jacobi preconditioning for A, maybe some row sum to form diagonal
// like in SIMPLEC to incorporate more of the off diagonal info (important in convection dominated flows)

class CVBC_AinvPC : public mfem::Solver
{
private:
   // CV bcs
   std::vector<ResistanceBC_T> *resistance_bcs = nullptr;
   bool resistance_bc = false;

   std::vector<WindkesselBC_T> *windkessel_bcs = nullptr;
   bool windkessel_bc = false;

   // diagonal of A
   mfem::Vector *A3d_diag = nullptr;
   mfem::Vector *Atotal_diag = nullptr;
   mfem::Vector *Abc_diag = nullptr;

   // velocity space
   mfem::FiniteElementSpace *vfes = nullptr;

public:
   CVBC_AinvPC(int truesize, mfem::Vector *A3d_diag_, mfem::Vector *Abc_diag_, mfem::Vector *Atotal_diag_, mfem::FiniteElementSpace *vfes_) : mfem::Solver(truesize), A3d_diag(A3d_diag_), Abc_diag(Abc_diag_), Atotal_diag(Atotal_diag_),
   vfes(vfes_) { };

   virtual void SetOperator(const mfem::Operator &op);

   void SetResistanceBC(std::vector<ResistanceBC_T> *resistance_bcs_)
   {
      resistance_bc = true;
      resistance_bcs = resistance_bcs_;
   }

   void SetWindkesselBC(std::vector<WindkesselBC_T> *windkessel_bcs_)
   {
      windkessel_bc = true;
      windkessel_bcs = windkessel_bcs_;
   }

   /// Operator application
   virtual void Mult (const mfem::Vector & x, mfem::Vector & y) const;

   /// Action of the transpose operator
   virtual void MultTranspose (const mfem::Vector & x, mfem::Vector & y) const;

   void GetTotalDiag(mfem::Vector &d)
   {
      d = *Atotal_diag;
   }
   // mfem::Vector &GetTotalDiag()
   // {
   //    return *Atotal_diag;
   // }

   mfem::Vector &GetAbcDiag()
   {
      return *Abc_diag;
   }

   virtual ~CVBC_AinvPC() {};
};