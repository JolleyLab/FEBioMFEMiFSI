#include "BodyForceIntegrator.hpp"

using namespace mfem;
using namespace std;
using namespace ifs_solver;

void BodyForceIntegrator::AssembleRHSElementVect(
   const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
{
   int vdim = Q.GetVDim();
   int dof  = el.GetDof();

   double val,cf;

   shape.SetSize(dof);       // vector of size dof

   elvect.SetSize(dof * vdim);
   elvect = 0.0;

   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // int intorder = 2*el.GetOrder() + Tr.OrderW();
      int intorder = 2*el.GetOrder();
      ir = &IntRules.Get(el.GetGeomType(), intorder);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      Tr.SetIntPoint (&ip);
      val = Tr.Weight();
      if (coeff)
      {
         val *= coeff->Eval(Tr,ip);
      }

      el.CalcShape(ip, shape);
      Q.Eval (Qvec, Tr, ip);

      for (int k = 0; k < vdim; k++)
      {
         cf = val * Qvec(k);

         for (int s = 0; s < dof; s++)
         {
            elvect(dof*k+s) += ip.weight * cf * shape(s);
         }
      }
   }
}