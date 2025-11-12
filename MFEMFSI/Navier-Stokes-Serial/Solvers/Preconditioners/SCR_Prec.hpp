#pragma once

// residual operator
#include "mfem.hpp"
#include "BlockPreconditioner.hpp"

using namespace std;

#ifdef MFEM_USE_MPI
// schur complement operator in matrix free fashion
class SchurOpSerial : public mfem::Operator
{
private:
   mfem::HypreParMatrix *C = nullptr; // not owned
   mfem::HypreParMatrix *D = nullptr; // not owned
   mfem::HypreParMatrix *G = nullptr; // not owned
   double Gcoeff, Dcoeff;
   mfem::Solver *Ainv = nullptr; // not owned

   bool equal_order = false;

   mutable mfem::Vector temp_y_2, temp_y;

public:
   SchurOpSerial(int size_p, int size_v, bool equal_order_) : mfem::Operator(size_p), equal_order(equal_order_)
   { 
      temp_y.SetSize(size_v); // square operator
      temp_y = 0.0;

      temp_y_2.SetSize(size_v); // square operator
      temp_y_2 = 0.0;
   };

   void SetAinvSolver(mfem::Solver *Ainv_)
   {
      Ainv = Ainv_;
   };

   void SetOperators(const mfem::Operator &D_, double &Dcoeff_, const mfem::Operator &G_, double &Gcoeff_, const mfem::Operator &C_)
   {
      // D  
      Dcoeff = Dcoeff_;
      const mfem::HypreParMatrix *new_D = dynamic_cast<const mfem::HypreParMatrix *>(&D_);
      MFEM_VERIFY(new_D, "new Operator must be a HypreParMatrix!");
      D = const_cast<mfem::HypreParMatrix *>(new_D);

      // G
      Gcoeff = Gcoeff_;
      const mfem::HypreParMatrix *new_G = dynamic_cast<const mfem::HypreParMatrix *>(&G_);
      MFEM_VERIFY(new_G, "new Operator must be a HypreParMatrix!");
      G = const_cast<mfem::HypreParMatrix *>(new_G);

      // C
      if (equal_order)
      {
         const mfem::HypreParMatrix *new_C = dynamic_cast<const mfem::HypreParMatrix *>(&C_);
         MFEM_VERIFY(new_C, "new Operator must be a HypreParMatrix!");
         C = const_cast<mfem::HypreParMatrix *>(new_C);
      }
   }

   virtual void Mult(const mfem::Vector &k, mfem::Vector &y) const override
   {
      // Apply S = C - D Ainv G

      // G->Mult(k,temp_y);
      G->Mult(Gcoeff,k,0.0,temp_y);

      Ainv->Mult(temp_y,temp_y_2);

      if (equal_order)
      {
         // D->Mult(temp_y_2,y);
         D->Mult(Dcoeff,temp_y_2,0.0,y);
         C->Mult(1.0,k,-1.0,y);
      }
      else
      {
         D->Mult(-Dcoeff,temp_y_2,0.0,y);
      }
   };

   virtual ~SchurOpSerial() noexcept {};
};

// Implements the SCR Method Preconditioner for block systems
//        Matrix System = [ A   G ]
//                        [ D   C ]
// where A = mass + diffusion + convection or just mass + diffusion
// G and D are the pressure gradient and velocity divergence terms
// C = 0 for inf-sup stable elements, C != 0 for equal order elements
//
//        Pscr =    [ I         0 ] [ A   0 ] [ I   A^-1 * G ]
//                  [ D*A^-1   I  ] [ 0   S ] [ 0        I   ]
//
// where S = C - DA^-1 * G, and its inverse
//
//   (Pscr)inv = [ I    -A^-1 * G ] [ A^-1   0 ] [ I        0 ]
//               [ 0         I    ] [ 0   S^-1 ] [ -DA^-1   I ]
//
//
class SCR_Prec_Serial : public BlockPreconditioner
{
private:
   // temp var to apply ops
   mutable mfem::Vector bhatp,yhatu;

   // matrix free schur complement
   SchurOpSerial *schur_comp = nullptr; // owned

   //Preconditioners
   mfem::Solver *A_inner_Inv = nullptr; // owned

   // options
   double a_inner_atol = 1E-10;
   double a_inner_rtol = 1E-6;

public:
   SCR_Prec_Serial(mfem::Array<int> &offsets_, mfem::Array<mfem::FiniteElementSpace *> &fes, bool equal_order_);

   void SetASolver(mfem::BlockOperator *jacobian) override;
   void SetSchurSolver(mfem::BlockOperator *jacobian) override;
   void FormSchurOperator(mfem::BlockOperator *jacobian) override;

   /// Operator application
   void Mult(const mfem::Vector & x, mfem::Vector & y) const override;

   void SetA_Inner_invTol(double a_i_atol_, double a_i_rtol_)
   {
      a_inner_atol = a_i_atol_;
      a_inner_rtol = a_i_rtol_;
   };

   virtual ~SCR_Prec_Serial();
};
#endif