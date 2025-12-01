#include "SIMPLE.hpp"
#include "../../../Utilities/sparsemat_ext.hpp"
#include "../../../Utilities/extended_operators.hpp"
#include "../../../Utilities/extended_solvers.hpp"
#include "../../../Utilities/hypre_func_extensions.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

SIMPLE_Prec::SIMPLE_Prec(Array<int> &offsets_,Array<FiniteElementSpace *> &fes, bool equal_order_) 
: BlockPreconditioner(offsets_,fes,equal_order_)
{
   tmp1.SetSize(offsets[1]-offsets[0]);
   tmp1 = 0.0;
   tmp2.SetSize(offsets[1]-offsets[0]);
   tmp2 = 0.0;
   tmp3.SetSize(offsets[2]-offsets[1]);
   tmp3 = 0.0;
   tmp4.SetSize(offsets[2]-offsets[1]);
   tmp4 = 0.0;
}

void SIMPLE_Prec::FormSchurOperator(mfem::BlockOperator *jacobian)
{
   // form Sparse approx of S = C - D(diagA)^-1 * G, formed each Newton iteration
   if (cv_bcs)
   {
      cvbc_PC->GetTotalDiag(Ad);
      // Ad = cvbc_PC->GetTotalDiag();
   }
   else
   {
      #ifdef MFEM_USE_MPI
      tempA = (HypreParMatrix *) &jacobian->GetBlock(0,0);
      #else
      tempA = (SparseMatrix *) &jacobian->GetBlock(0,0);
      #endif
      tempA->GetDiag(Ad);
   }   
   
   if(tempS) delete tempS;
   // any coeffs on G and D blocks, all other blocks will have correct factors already
   Gcoeff = jacobian->GetBlockCoef(0,1);
   Dcoeff = jacobian->GetBlockCoef(1,0);

   #ifdef MFEM_USE_MPI
   tempG = (HypreParMatrix *) &jacobian->GetBlock(0,1);   
   HypreParMatrixInvScaleRows(tempG,Ad,Gcoeff*Dcoeff*-1.0); // G = G_0 * 1/diag(A) * (af * af * -1.0)
   #else
   tempG = (SparseMatrix *) &jacobian->GetBlock(0,1);
   SparseMatInvScaleRows(*tempG,Ad,Gcoeff*Dcoeff*-1.0); // G = G_0 * 1/diag(A) * (af * af * -1.0)
   #endif
   
   #ifdef MFEM_USE_MPI
   tempD = (HypreParMatrix *) &jacobian->GetBlock(1,0);
   tempS = ParMult(tempD,tempG,true);
   HypreParMatrixScaleRows(tempG,Ad,1.0/(Gcoeff*Dcoeff*-1.0)); // G = [G_0 * 1/diag(A) * (af * af * -1.0)] * [diag(A) / (af * af * -1.0)] = G_0
   #else
   tempD = (SparseMatrix *) &jacobian->GetBlock(1,0);
   tempS = mfem::Mult(*tempD,*tempG); // tempS is owned by this class
   SparseMatScaleRows(*tempG,Ad,1.0/(Gcoeff*Dcoeff*-1.0)); // G = [G_0 * 1/diag(A) * (af * af * -1.0)] * [diag(A) / (af * af * -1.0)] = G_0
   #endif

   if (equal_order)
   {
      if (S) delete S;
      
      #ifdef MFEM_USE_MPI
      tempC = (HypreParMatrix *) &jacobian->GetBlock(1,1);
      S = ParAdd(tempS,tempC);
      #else
      tempC = (SparseMatrix *) &jacobian->GetBlock(1,1);
      S = Add(*tempS,*tempC); // S is owned by this class
      #endif

      SInv->SetOperator(*S);      
   }
   else
   {
      SInv->SetOperator(*tempS); 
   }

   // modify Sinv operator to account for bcs (elim any zero rows)
   if (equal_order)
   {
      S->EliminateZeroRows();
   }
   else
   {
      tempS->EliminateZeroRows();
   }
}

// Operator application
void SIMPLE_Prec::Mult (const Vector & x,Vector & y) const
{   
   MFEM_ASSERT(x.Size() == width, "incorrect input Vector size");
   MFEM_ASSERT(y.Size() == height, "incorrect output Vector size");

   yblock.Update(y.GetData(),offsets); // [u p]
   xblock.Update(x.GetData(),offsets); // [bu bp]
   // y = 0.0;

   // Step 1
   AInv->Mult(xblock.GetBlock(0),tmp1);
   #ifdef MFEM_USE_MPI
   tempD->Mult(-Dcoeff,tmp1,0.0,tmp3);
   #else
   tmp3 = 0.0;
   tempD->AddMult(tmp1,tmp3,-Dcoeff);
   #endif
   tmp3 += xblock.GetBlock(1);

   // Step 2
   // tmp3 = -D*Ainv*bu + bp
   // tmp1 = Ainv*bu
   SInv->Mult(tmp3,tmp4);

   // Step 3
   // tmp4 = Sinv(-D*Ainv*bu + bp)
   // tmp1 = Ainv*bu
   // tempG = (diagA)^-1 * G
   #ifdef MFEM_USE_MPI
   tempG->Mult(Gcoeff,tmp4,0.0,tmp2);
   #else
   tmp2 = 0.0;
   tempG->AddMult(tmp4,tmp2,Gcoeff);
   #endif
   tmp2 /= Ad;
   add(1.0,tmp1,-1.0,tmp2,yblock.GetBlock(0));
   yblock.GetBlock(1).Set(1.0/alpha,tmp4);
}

SIMPLE_Prec::~SIMPLE_Prec()
{
   
}