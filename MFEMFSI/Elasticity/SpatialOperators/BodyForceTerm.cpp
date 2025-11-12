#include "BodyForceTerm.hpp"

#include <iomanip>

using namespace mfem;
using namespace std;

void VectorDomainLFIntegrator_Ref::AssembleRHSElementVect(
   const FiniteElement &el, ElementTransformation &Tr, Vector &elvect)
{
   int vdim = Q.GetVDim();
   int dof  = el.GetDof();

   shape.SetSize(dof);
   vec.SetSize(vdim);

   elvect.SetSize(dof * vdim);
   elvect = 0.0;

   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      int intorder = 2*el.GetOrder() + Tr.OrderW();
      ir = &IntRules.Get(el.GetGeomType(), intorder);
   }

   double ws = 0.0;

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);

      // build transformation if necessary
      if (formulation_type == form_type::UpdatedLagrangian)
      {
         // always compute integral in reference config
         BuildInitTransformation(Tr,ip,&el);
         Q.Eval(vec, *init_trans, ip);
         vec *= init_trans->Weight() * ip.weight;

         if (rho)
         {
            vec *= rho->Eval(*init_trans,ip);
         }
      }
      else if (formulation_type == form_type::TotalLagrangian)
      {
         Q.Eval(vec, Tr, ip);
         vec *= Tr.Weight() * ip.weight;

         if (rho)
         {
            vec *= rho->Eval(Tr,ip);
         }
      }
      else
      {
         MFEM_ABORT("Error: unsupported formulation type");
      }
      
      el.CalcShape(ip, shape);
      for (int k = 0; k < vdim; k++)
      {
         for (int s = 0; s < dof; s++)
         {
            elvect(dof*k+s) += vec(k) * shape(s);
         }
      }
   }
}

void VectorDomainLFIntegrator_Ref::AssembleRHSElementVect(
   const FiniteElement &el, FaceElementTransformations &Tr, Vector &elvect)
{
   int vdim = Q.GetVDim();
   int dof  = el.GetDof();

   shape.SetSize(dof);
   vec.SetSize(vdim);

   elvect.SetSize(dof * vdim);
   elvect = 0.0;

   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int intorder = 2*el.GetOrder();
      ir = &IntRules.Get(Tr.GetGeometryType(), intorder);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      // Use Tr transformation in case Q depends on boundary attribute
      // build transformation if necessary
      if (formulation_type == form_type::UpdatedLagrangian)
      {
         // always compute integral in reference config
         BuildInitTransformation(Tr,ip,&el);
         Q.Eval(vec, *init_trans, ip);
         vec *= init_trans->Weight() * ip.weight;

         if (rho)
         {
            vec *= rho->Eval(*init_trans,ip);
         }
      }
      else if (formulation_type == form_type::TotalLagrangian)
      {
         Q.Eval(vec, Tr, ip);
         vec *= Tr.Weight() * ip.weight;

         if (rho)
         {
            vec *= rho->Eval(Tr,ip);
         }
      }
      else
      {
         MFEM_ABORT("Error: unsupported formulation type");
      }      

      el.CalcShape(eip, shape);
      for (int k = 0; k < vdim; k++)
      {
         for (int s = 0; s < dof; s++)
         {
            elvect(dof*k+s) += vec(k) * shape(s);
         }
      }
   }
}

// void VectorBoundaryLFIntegrator_Ref::BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip, const FiniteElement *curr_elem)
// {
//    // setup transformation from parent element -> initial configuration
//    // only rebuild if T.ElementNo != init_trans.ElementNo
//    if (init_trans.ElementNo != T.ElementNo)
//    {
//       init_trans.SetFE(curr_elem);

//       // get original node positions
//       int dim = curr_elem->GetDim();
//       Vector vals, x_orig;
//       Array<int> vdofs;

      // // build point matrix for init_trans
      // DenseMatrix init_pm;
      // int idx, idx1, idx2;
      // int vdim_num;      

      // // get vdofs depending on type of transformation
      // if (T.ElementType == ElementTransformation::BDR_ELEMENT)
      // {
      //    // boundary element transformation
      //    auto dof_trans = fespace->GetBdrElementVDofs(T.ElementNo,vdofs);

      //    // set parameters
      //    vdim_num = dim + 1;
      //    x_orig.SetSize(vdim_num);
      //    init_pm.SetSize(vdim_num,curr_elem->GetDof());
      // }
      // else if (T.ElementType == ElementTransformation::ELEMENT)
      // {
      //    // element transformation
      //    auto dof_trans = fespace->GetElementVDofs(T.ElementNo,vdofs);

      //    // set parameters
      //    vdim_num = dim;
      //    x_orig.SetSize(vdim_num);
      //    init_pm.SetSize(vdim_num,curr_elem->GetDof());
      // }
      // else
      // {
      //    MFEM_ABORT("TO DO: ERROR, rebuilding init_bdr_transformation, not implemented element type");
      // }
      
//       orig_nodes->GetSubVector(vdofs,vals);
      // int shift = vdofs.Size()/vdim_num;
      // for (int j=0; j<vdim_num; j++)
      // {
      //    idx = j;
      //    if (vdim_num == 2) idx1 = j + shift;
      //    if (vdim_num == 3) idx2 = j + 2*shift;

      //    x_orig(0) = vals(idx);
      //    if (vdim_num == 2) x_orig(1) = vals(idx1);
      //    if (vdim_num == 3) x_orig(2) = vals(idx2);

      //    // set colum to x_orig
      //    init_pm.SetCol(j,x_orig);
      // }
//       // init_pm.Print();
//       // copy over T details
//       init_trans.SetPointMat(init_pm);
//       init_trans.ElementNo = T.ElementNo;
//       init_trans.ElementType = T.ElementType;
//    }
//    // always set current int point
//    init_trans.SetIntPoint(&ip);
// }

void VectorDomainLFIntegrator_Ref::BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip, const FiniteElement *curr_elem)
{
   auto curr_trans = dynamic_cast<IsoparametricTransformation *>(&T);
   // auto curr_elem = curr_trans->GetFE();

   // setup transformation from parent element -> initial configuration
   // only rebuild if T.ElementNo != init_trans.ElementNo
   if (init_trans == nullptr || init_trans->ElementNo != T.ElementNo)
   {
      if (init_trans) delete init_trans;
      init_trans = new IsoparametricTransformation();
      // if(!init_trans) init_trans = new IsoparametricTransformation();
      init_trans->SetFE(curr_elem);

      // get original node positions
      int dim = curr_elem->GetDim();
      Vector vals, x_orig;
      x_orig = 0.0;
      Array<int> vdofs;
      
      // build point matrix for init_trans
      DenseMatrix init_pm;
      int idx, idx1, idx2;
      int vdim_num;      

      // get vdofs depending on type of transformation
      if (T.ElementType == ElementTransformation::BDR_ELEMENT)
      {
         // boundary element transformation
         auto dof_trans = fespace->GetBdrElementVDofs(T.ElementNo,vdofs);

         // set parameters
         vdim_num = dim + 1;
         x_orig.SetSize(vdim_num);
         init_pm.SetSize(vdim_num,curr_elem->GetDof());
      }
      else if (T.ElementType == ElementTransformation::ELEMENT)
      {
         // element transformation
         auto dof_trans = fespace->GetElementVDofs(T.ElementNo,vdofs);

         // set parameters
         vdim_num = dim;
         x_orig.SetSize(vdim_num);
         init_pm.SetSize(vdim_num,curr_elem->GetDof());
      }
      else
      {
         MFEM_ABORT("TO DO: ERROR, rebuilding init_bdr_transformation, not implemented element type");
      }      
      orig_nodes->GetSubVector(vdofs,vals);
      int nd = vdofs.Size()/vdim_num;
      for (int j=0; j<nd; j++)
      {
         idx = j;
         idx1 = j + nd;
         if (vdim_num == 3) idx2 = j + 2*nd;

         x_orig(0) = vals(idx);
         x_orig(1) = vals(idx1);
         if (vdim_num == 3) x_orig(2) = vals(idx2);

         // set colum to x_orig
         init_pm.SetCol(j,x_orig);
      }
      // init_pm.Print();
      // copy over T details
      init_trans->SetPointMat(init_pm);
      init_trans->ElementNo = T.ElementNo;
      init_trans->ElementType = T.ElementType;
   }
   // always set current int point
   init_trans->SetIntPoint(&ip);
}