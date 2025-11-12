#pragma once

#include "mfem.hpp"

using namespace mfem;
using namespace std;

// block preconditioner for incompressible elasticity, from ex19p on MFEM
// Custom block preconditioner for the Jacobian of the incompressible nonlinear
// elasticity operator. It has the form
//
// P^-1 = [ K^-1 0 ][ I -B^T ][ I  0           ]
//        [ 0    I ][ 0  I   ][ 0 -\gamma S^-1 ]
//
// where the original Jacobian has the form
//
// J = [ K B^T ]
//     [ B 0   ]
//
// and K^-1 is an approximation of the inverse of the displacement part of the
// Jacobian and S^-1 is an approximation of the inverse of the Schur
// complement S = B K^-1 B^T. The Schur complement is approximated using
// a mass matrix of the pressure variables.
class JacobianPreconditioner : public Solver
{
protected:
   // Finite element spaces for setting up preconditioner blocks
   Array<ParFiniteElementSpace *> spaces;

   // Offsets for extracting block vector segments
   Array<int> &block_trueOffsets;

   // Jacobian for block access
   BlockOperator *jacobian;

   // Scaling factor for the pressure mass matrix in the block preconditioner
   double gamma;

   // Objects for the block preconditioner application
   ParBilinearForm *PressureMassForm = nullptr;
   HypreParMatrix *pressure_mass = nullptr;
   Solver *mass_solver = nullptr;
   Solver *mass_prec = nullptr;
   Solver *stiff_solver = nullptr;
   Solver *stiff_prec = nullptr;

public:
   JacobianPreconditioner(Array<ParFiniteElementSpace *> &fes, Array<int> &offsets);

   virtual void Mult(const Vector &k, Vector &y) const;
   virtual void SetOperator(const Operator &op);

   virtual ~JacobianPreconditioner();
};