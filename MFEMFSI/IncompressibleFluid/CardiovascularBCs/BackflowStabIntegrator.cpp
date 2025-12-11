#include "BackflowStabIntegrator.hpp"

using namespace mfem;
using namespace std;

void BackflowStabIntegrator::AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &T,
                                      const Vector &elfun,
                                      Vector &elvect)
{
   // note: elfun is list of dofs, or fem coefficients

   // here element transformation and element is on boundary
   int dim = el.GetDim()+1;
   int dof = el.GetDof();

   shape.SetSize(dof);
   nor.SetSize(dim);
   u.SetSize(dim);
   un_u.SetSize(dim);
   u_bf.SetSize(dim);
   if (ale) uhat.SetSize(dim);
   EF.UseExternalData(elfun.GetData(),dof,dim);
   elvect.SetSize(dof * dim);
   elvect = 0.0;

   double w = 0.0;

   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      int order = 2*el.GetOrder();
      ir = &IntRules.Get(el.GetGeomType(), order);
   }
   
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      // set the current quad point in reference space
      const IntegrationPoint &ip = ir->IntPoint(i);
      T.SetIntPoint(&ip);

      // compute normal
      if (dim > 1)
      {
         CalcOrtho(T.Jacobian(), nor);
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

      // evaluate basis in reference space
      el.CalcShape(ip,shape);
      // compute u
      EF.MultTranspose(shape,u);

      // compute (u - uhat) * n
      u_bf = u;
      if (ale)
      {
         vec_coeff->Eval(uhat,T,ip);
         u_bf -= uhat;
      }
      double u_dot_n = u_bf * nor;
      // check if backflow
      // u_dot_n = u_dot_n < 0.0 ? u_dot_n : 0.0;
      u_dot_n = std::min(u_dot_n, 0.0);

      // (u - uhat) dot n * u
      un_u.Set(u_dot_n,u);

      w = ip.weight * T.Weight();
      if (coeff) { w *= coeff->Eval(T,ip); }
      w *= beta;

      for (int a = 0; a < dof; a++)
      {
         elvect(a) += w * shape(a) * un_u(0);
         elvect(a + dof) += w * shape(a) * un_u(1);
         if (dim == 3) elvect(a + 2*dof) += w * shape(a) * un_u(2);
      }
   }
}

void BackflowStabIntegrator::AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &trans,
                                     const Vector &elfun,
                                     DenseMatrix &elmat)
{
   // here element transformation and element is on boundary
   int dim = el.GetDim()+1;
   int dof = el.GetDof();

   shape.SetSize(dof);
   nor.SetSize(dim);
   u.SetSize(dim);
   un_u.SetSize(dim);
   u_bf.SetSize(dim);
   if (ale) uhat.SetSize(dim);
   EF.UseExternalData(elfun.GetData(),dof,dim);
   
   elmat.SetSize((dof)*dim);
   elmat = 0.0;

   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      int order = 2*el.GetOrder();
      ir = &IntRules.Get(el.GetGeomType(), order);
   }
   
   double w = 0.0;
   a11.SetSize(dof);
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      // set the current quad point in reference space
      const IntegrationPoint &ip = ir->IntPoint(i);
      trans.SetIntPoint(&ip);

      // compute normal
      if (dim > 1)
      {
         CalcOrtho(trans.Jacobian(), nor);
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

      // evaluate basis in reference space
      el.CalcShape(ip,shape);
      // compute u
      EF.MultTranspose(shape,u);

      // compute (u - uhat) * n
      u_bf = u;
      if (ale)
      {
         vec_coeff->Eval(uhat,trans,ip);
         u_bf -= uhat;
      }
      double u_dot_n = u_bf * nor;
      // check if backflow
      // u_dot_n = u_dot_n < 0.0 ? u_dot_n : 0.0;
      u_dot_n = std::min(u_dot_n, 0.0);

      w = ip.weight*trans.Weight();
      if (coeff) { w *= coeff->Eval(trans,ip); } // for density
      w *= u_dot_n * beta;

      for (int a=0; a<dof; ++a)
      {
         for (int b=0; b<dof; ++b)
         {
            double bf_val = w * shape(a) * shape(b);
            elmat(a,b)                               += bf_val;
            elmat(a + dof,b + dof)                   += bf_val;
            if (dim == 3) elmat(a + 2*dof,b + 2*dof) += bf_val;
         }
      }
   }
}

/////////////////////////////////
// block form

void BlockBackflowStabIntegrator::AssembleElementVector(const Array<const FiniteElement *> &el,
                                       ElementTransformation &Tr,
                                       const Array<const Vector *> &elfun,
                                       const Array<Vector *> &elvec)
{
   // here element transformation and element is on boundary
   int dim = el[0]->GetDim()+1;
   int dof = el[0]->GetDof();

   shape.SetSize(dof);
   nor.SetSize(dim);
   u.SetSize(dim);
   un_u.SetSize(dim);
   u_bf.SetSize(dim);
   if (mesh_velo) uhat.SetSize(dim);
   EF.UseExternalData(elfun[0]->GetData(),dof,dim);
   
   elvec[0]->SetSize(dof*dim);
   (*elvec[0]) = 0.0;
   
   elvec[1]->SetSize(0);
   // only size and use the first block, the second will be ignored!
   double w = 0.0;

   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      // int order = el[0]->GetOrder()*3 + Tr.OrderW();
      // ir = &IntRules.Get(el[0]->GetGeomType(), order);
      // int order = 2*el[0]->GetOrder() + 3;
      int order = 2*el[0]->GetOrder(); // exact integration for linear terms
      ir = &IntRules.Get(el[0]->GetGeomType(), order);
   }
   
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      // set the current quad point in reference space
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);

      // compute normal
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

      // evaluate basis in reference space
      el[0]->CalcShape(ip,shape);
      // compute u
      EF.MultTranspose(shape,u);

      // compute (u - uhat) * n
      u_bf = u;
      if (ale)
      {
         // mesh_velo->GetVectorValue(Tr,ip,uhat);
         mesh_velo->Eval(uhat,Tr,ip);
         u_bf -= uhat;
      }
      double u_dot_n = u_bf * nor;
      // check if backflow
      // u_dot_n = u_dot_n < 0.0 ? u_dot_n : 0.0;
      u_dot_n = std::min(u_dot_n, 0.0);
      // if (u_dot_n == 0.0) continue;      

      // (u - uhat) dot n * u
      un_u.Set(u_dot_n,u);

      w = ip.weight * Tr.Weight();
      if (coeff) { w *= coeff->Eval(Tr,ip); }
      w *= -beta;

      for (int a = 0; a < dof; a++)
      {
         (*elvec[0])(a) += w * shape(a) * un_u(0);
         (*elvec[0])(a + dof) += w * shape(a) * un_u(1);
         if (dim == 3) (*elvec[0])(a + 2*dof) += w * shape(a) * un_u(2);
      }
   }

   // cout << "elvec = " ;
   // (*elvec[0]).Print(cout,10);
}

void BlockBackflowStabIntegrator::AssembleElementGrad(const Array<const FiniteElement*> &el,
                                     ElementTransformation &Tr,
                                     const Array<const Vector *> &elfun,
                                     const Array2D<DenseMatrix *> &elmats)
{
   // here element transformation and element is on boundary
   int dim = el[0]->GetDim()+1;
   int dof = el[0]->GetDof();

   shape.SetSize(dof);
   nor.SetSize(dim);
   u.SetSize(dim);
   un_u.SetSize(dim);
   u_bf.SetSize(dim);
   if (mesh_velo) uhat.SetSize(dim);
   EF.UseExternalData(elfun[0]->GetData(),dof,dim);

   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      // int order = el[0]->GetOrder()*3 + Tr.OrderW();
      // ir = &IntRules.Get(el[0]->GetGeomType(), order);
      // int order = 2*el[0]->GetOrder() + 3;
      int order = 2*el[0]->GetOrder(); // exact integration for linear terms
      ir = &IntRules.Get(el[0]->GetGeomType(), order);
   }

   double w = 0.0;
   elmats(0,0)->SetSize(dof*dim, dof*dim); // only size (0,0) block, others will be ignored!
   *elmats(0,0) = 0.0;

   elmats(0,1)->SetSize(0);
   elmats(1,0)->SetSize(0);
   elmats(1,1)->SetSize(0);

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      // set the current quad point in reference space
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);

      // compute normal
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

      // evaluate basis in reference space
      el[0]->CalcShape(ip,shape);
      // compute u
      EF.MultTranspose(shape,u);

      // compute (u - uhat) * n
      u_bf = u;
      if (ale)
      {
         // mesh_velo->GetVectorValue(Tr,ip,uhat);
         mesh_velo->Eval(uhat,Tr,ip);
         u_bf -= uhat;
      }
      double u_dot_n = u_bf * nor;
      // check if backflow
      // u_dot_n = u_dot_n < 0.0 ? u_dot_n : 0.0;
      u_dot_n = std::min(u_dot_n, 0.0);    
      // if (u_dot_n == 0.0) continue; 

      w = ip.weight*Tr.Weight();
      if (coeff) { w *= coeff->Eval(Tr,ip); } // for density
      w *= u_dot_n * -1.0 * beta;
      double af = gen_alpha ? alpha_f*gamma*(*dt) : 1.0;
      w *= af;

      // 11 block
      for (int a=0; a<dof; ++a)
      {
         for (int b=0; b<dof; ++b)
         {
            double bf_val = w * shape(a) * shape(b);
            (*elmats(0,0))(a,b)                               += bf_val;
            (*elmats(0,0))(a + dof,b + dof)                   += bf_val;
            if (dim == 3) (*elmats(0,0))(a + 2*dof,b + 2*dof) += bf_val;
         }
      }
   }
}