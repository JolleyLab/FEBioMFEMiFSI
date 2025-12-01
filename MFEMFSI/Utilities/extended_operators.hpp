#pragma once

#include "mfem.hpp"
#include "../Navier-Stokes/CardiovascularBCs/CV_Op.hpp"

using namespace ins_solver;

// sum operator for tangent stiffness matrix A block with CV bcs
class SumOp : public mfem::Operator
{
private:
   const mfem::Operator *A = nullptr; // not owned (A3d)
   const CVBC_Op *B = nullptr; // not owned (cv)
   const mfem::Operator *C = nullptr; // not owned (stabilization)

   const mfem::Array<int> *vel_ess_tdofs = nullptr; // not owned

   double cv_time_factor;

   mutable mfem::Vector temp;
public:
   SumOp(int size, const mfem::Array<int> &vel_ess_tdofs_) : mfem::Operator(size), vel_ess_tdofs(&vel_ess_tdofs_) 
   {
      temp.SetSize(size); // square operators!
      temp = 0.0;
   };

   void SetOps(const mfem::Operator &A_, const CVBC_Op &B_)
   {
      A = &A_;
      B = &B_;
   }
   
   void SetCV_timefactor(double &af)
   {
      cv_time_factor = af;
   }

   virtual void Mult(const mfem::Vector &k, mfem::Vector &y) const override;

   void SetC(const mfem::Operator &C_)
   {
      C = &C_;
   }

   const mfem::Operator *GetA() const
   {
      return A;
   }

   virtual ~SumOp() {};
};