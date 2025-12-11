#include "aSCR_Prec.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

#ifdef MFEM_USE_MPI
aSCR_Prec::aSCR_Prec(Array<int> &offsets_,Array<FiniteElementSpace *> &fes, bool equal_order_) 
: BlockPreconditioner(offsets_,fes,equal_order_)
{
   bhatp.SetSize(offsets[2]-offsets[1]);
   bhatp = 0.0;
   yhatu.SetSize(offsets[1]-offsets[0]);
   yhatu = 0.0;

   tmp1.SetSize(offsets[2]-offsets[1]);
   tmp1 = 0.0;
   tmp2.SetSize(offsets[1]-offsets[0]);
   tmp2 = 0.0;
   tmp3.SetSize(offsets[1]-offsets[0]);
   tmp3 = 0.0;
}

void aSCR_Prec::FormSchurOperator(mfem::BlockOperator *jacobian)
{  
   // form Sparse approx of S = C - D(diagA)^-1 * G, formed each Newton iteration
   if (cv_bcs)
   {
      cvbc_PC->GetTotalDiag(Ad);
   }
   else
   {
      tempA = (HypreParMatrix *) &jacobian->GetBlock(0,0);
      tempA->GetDiag(Ad);
   }

   if(tempS) delete tempS;
   // any coeffs on G and D blocks, all other blocks will have correct factors already
   Gcoeff = jacobian->GetBlockCoef(0,1);
   Dcoeff = jacobian->GetBlockCoef(1,0);

   tempG = (HypreParMatrix *) &jacobian->GetBlock(0,1);
   HypreParMatrixInvScaleRows(tempG,Ad,Gcoeff*Dcoeff*-1.0); // G = G_0 * 1/diag(A) * (af * af * -1.0)

   tempD = (HypreParMatrix *) &jacobian->GetBlock(1,0);
   tempS = ParMult(tempD,tempG,true);

   HypreParMatrixScaleRows(tempG,Ad,1.0/(Gcoeff*Dcoeff*-1.0)); // G = [G_0 * 1/diag(A) * (af * af * -1.0)] * [diag(A) / (af * af * -1.0)] = G_0

   if (equal_order)
   {
      if (S) delete S;
      
      tempC = (HypreParMatrix *) &jacobian->GetBlock(1,1);
      S = ParAdd(tempS,tempC);
      SInv->SetOperator(*S);      
   }
   else
   {
      SInv->SetOperator(*tempS); 
   }

   // static_cast<HypreParMatrix *>(tempS)->EliminateZeroRows(); // this is a little bit faster than above!
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
void aSCR_Prec::Mult(const Vector & x,Vector & y) const
{
   MFEM_ASSERT(x.Size() == width, "incorrect input Vector size");
   MFEM_ASSERT(y.Size() == height, "incorrect output Vector size");

   yblock.Update(y.GetData(),offsets); // [u p]
   xblock.Update(x.GetData(),offsets); // [bu bp]

   // Step 1
   AInv->Mult(xblock.GetBlock(0),yhatu);
   tempD->Mult(Dcoeff,yhatu,0.0,tmp1);
   // add(1.0,xblock.GetBlock(1),-1.0,tmp1,bhatp);
   add(1.0/alpha,xblock.GetBlock(1),-1.0/alpha,tmp1,bhatp);

   // Step 2
   SInv->Mult(bhatp,yblock.GetBlock(1));
   // yblock.GetBlock(1) *= 1.0/alpha;

   // Step 3
   tempG->Mult(Gcoeff,yblock.GetBlock(1),0.0,tmp2);
   AInv->Mult(tmp2,tmp3);
   add(1.0,yhatu,-1.0,tmp3,yblock.GetBlock(0));
}

aSCR_Prec::~aSCR_Prec()
{

}
#endif