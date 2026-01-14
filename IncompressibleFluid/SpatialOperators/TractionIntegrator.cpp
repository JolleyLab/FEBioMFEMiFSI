#include "TractionIntegrator.hpp"

using namespace mfem;
using namespace std;
using namespace ifs_solver;

void TractionIntegrator::AssembleRHSElementVect(
   const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
{
   int dim = el.GetDim() + 1; // see MFEM BoundaryNormalLFIntegrator!
   int dof = el.GetDof();
   Vector nor(dim), Qvec(dim);
   nor = 0.0;
   Qvec = 0.0;
   DenseMatrix gradu(dim);
   gradu = 0.0;

   shape.SetSize(dof);
   elvect.SetSize(dof*dim);
   elvect = 0.0;

   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // int order = 2*el.GetOrder() + Tr.OrderW();
      int order = 2*el.GetOrder(); // exact for linear terms
      ir = &IntRules.Get(Tr.GetGeometryType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetIntPoint(&ip);

      // normal, evaluate only for matrix coefficients or scalar pressures
      if (MQ || scalar_coeff)
      {
         if (dim > 1)
         {
            CalcOrtho(Tr.Jacobian(), nor);
         }
         else
         {
            nor[0] = 1.0;
         }

         // normalize nor, see example 18 in general it is not a unit normal
         double normag = 0;
         for (int ii = 0; ii < nor.Size(); ii++)
         {
            normag += nor(ii)*nor(ii);
         }
         normag = sqrt(normag);
         nor *= 1/normag;
      }
      
      // evaulate coeff
      if(Q)
      {
         Q->Eval(Qvec, Tr, ip);
      }
      else if(MQ)
      {
         MQ->Eval(gradu,Tr,ip);
         gradu.Mult(nor,Qvec);
      }
      else if (scalar_coeff)
      {
         double val = scalar_coeff->Eval(Tr,ip);
         if (pressure) val *= -1.0; // traction is -p*n
         Qvec.Set(val,nor);
      }
      else
      {
         MFEM_ABORT("internal error");
      }

      el.CalcShape(ip, shape);

      double w = ip.weight*Tr.Weight();
      if (coeff) { w *= coeff->Eval(Tr,ip); }
      
      for (int j = 0; j < dof; j++)
      {
         elvect(j) += w*Qvec(0)*shape(j);
         elvect(j + dof) += w*Qvec(1)*shape(j);
         if (dim == 3) elvect(j + 2*dof) += w*Qvec(2)*shape(j);
      }
   }
}