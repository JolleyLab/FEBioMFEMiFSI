#pragma once

// residual operator
#include "mfem.hpp"
#include "BlockPreconditioner.hpp"

// Implements the SIMPLE Preconditioner for block systems
//        Matrix System = [ A   G ]
//                        [ D   C ]
// where A = mass + diffusion + convection or just mass + diffusion
// G and D are the pressure gradient and velocity divergence terms
// C = 0 for inf-sup stable elements, C != 0 for equal order elements
//
//        Psimple = [ I         0 ] [ A   0 ] [ I   (diagA)^-1 * G ]
//                  [ D*A^-1   I  ] [ 0   S ] [ 0        alpha I   ]
//
// where S = C - D(diagA)^-1 * G, alpha \in (0,1], and its inverse
//
//   (Psimple)inv = [ I    -(diagA)^-1 * G ] [ A^-1   0 ] [ I        0 ]
//                  [ 0       1/alpha I    ] [ 0   S^-1 ] [ -DA^-1   I ]
//
//
class SIMPLE_Prec_Serial : public BlockPreconditioner
{
private:
   // options
   // double alpha = 0.5; // under relaxation parameter for pressure update
   double alpha = 0.9; // under relaxation parameter for pressure update
   // double alpha = 1.0; // no relaxation for pressure update // can cause issues sometimes

   mutable mfem::Vector tmp4;

public:
   SIMPLE_Prec_Serial(mfem::Array<int> &offsets_, mfem::Array<mfem::FiniteElementSpace *> &fes, bool equal_order_);

   void FormSchurOperator(mfem::BlockOperator *jacobian) override;

   void Mult(const mfem::Vector & x, mfem::Vector & y) const override;

   virtual ~SIMPLE_Prec_Serial();
};