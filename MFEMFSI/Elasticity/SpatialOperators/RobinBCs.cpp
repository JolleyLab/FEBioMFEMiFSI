#include "RobinBCs.hpp"

using namespace mfem;
using namespace std;

void RobinBCIntegrator::AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &T,
                                      const Vector &elfun,
                                      Vector &elvect)
{
   // here element transformation and element is on boundary
   int dim = el.GetDim()+1;
   int nd = el.GetDof();

   shape.SetSize(nd);
   u.SetSize(dim);
   udot.SetSize(dim);
   EF.UseExternalData(elfun.GetData(),nd,dim);
   
   elvect.SetSize(nd * dim);
   elvect = 0.0;

   double w = 0.0;

   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      int order = el.GetOrder()*2 + 3;
      ir = &IntRules.Get(el.GetGeomType(), order);
   }
   
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      // set the current quad point in reference space
      const IntegrationPoint &ip = ir->IntPoint(i);
      T.SetIntPoint(&ip);

      // compute shape
      el.CalcShape(ip,shape);
      // compute u
      EF.MultTranspose(shape,u);
      // compute udot
      curr_v->GetVectorValue(T,ip,udot);

      // compute -ks*u - cs*udot
      u *= -ks;
      u.Add(-cs,udot);

      w = ip.weight * T.Weight();
      for (int a = 0; a < nd; a++)
      {
         elvect(a) += w * shape(a) * u(0);
         elvect(a + nd) += w * shape(a) * u(1);
         if (dim == 3) elvect(a + 2*nd) += w * shape(a) * u(2);
      }
   }
}

void RobinBCIntegrator::AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &T,
                                     const Vector &elfun,
                                     DenseMatrix &elmat)
{
   // here element transformation and element is on boundary
   int dim = el.GetDim()+1;
   int nd = el.GetDof();

   shape.SetSize(nd);
   u.SetSize(dim);
   udot.SetSize(dim);
   EF.UseExternalData(elfun.GetData(),nd,dim);
   
   elmat.SetSize(nd * dim);
   elmat = 0.0;

   double w = 0.0;
   double val = 0.0;
   double af = alpha_f*beta*(*curr_dt)*(*curr_dt);
   double afv = alpha_f*gamma*(*curr_dt);

   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      int order = el.GetOrder()*2 + 3;
      ir = &IntRules.Get(el.GetGeomType(), order);
   }
   
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      // set the current quad point in reference space
      const IntegrationPoint &ip = ir->IntPoint(i);
      T.SetIntPoint(&ip);

      // compute shape
      el.CalcShape(ip,shape);

      w = ip.weight * T.Weight();
      if (div_factor) w /= af;
      for (int a = 0; a < nd; a++)
      {
         for (int b = 0; b < nd; b++)
         {
            val = shape(a)*shape(b)*w;

            elmat(a,b) += val*(-ks*af - cs*afv);
            elmat(a+nd,b+nd) += val*(-ks*af - cs*afv);
            if (dim == 3) elmat(a+2*nd,b+2*nd) += val*(-ks*af - cs*afv);
         }
      }
   }
}