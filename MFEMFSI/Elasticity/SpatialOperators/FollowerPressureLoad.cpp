#include "FollowerPressureLoad.hpp"

#include <iomanip>

using namespace mfem;
using namespace std;

// CG TERMS
void FollowerPressureLoad::AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &trans,
                                      const Vector &elfun,
                                      Vector &elvect)
{
   const int nd = el.GetDof();
   const int dim = el.GetDim() + 1; // for boundary elements

   Vector shape, nor(dim), temp_vec(nd);

   shape.SetSize(nd);
   shape = 0.0;
   nor = 0.0;
   temp_vec = 0.0;

   elvect.SetSize(nd*dim);
   elvect = 0.0;

   if (formulation_type == form_type::TotalLagrangian)
   {
      EF.UseExternalData(elfun.GetData(), nd, dim);
   }

   double w = 0.0;
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order = 2 * el.GetOrder() + trans.OrderW(); // correct order? Same as vectorboundaryintegrator
      ir = &IntRules.Get(trans.GetGeometryType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {     
      const IntegrationPoint &ip = ir->IntPoint(i);

      if (formulation_type == form_type::TotalLagrangian)
      {
         BuildCurrentTransformation(trans,ip,&EF);
         // curr_trans->SetIntPoint(&ip); // this is done inside above function
         
         w = ip.weight; // No detJ, cancels out for this integrator!
         w *= pcoeff->Eval(*curr_trans, ip); // p*w, evaluate thru face transformation
         if(coeff)
         {
            w *= coeff->Eval(*curr_trans,ip);
         }      

         // compute the normal (unnormalized!)
         CalcOrtho(curr_trans->Jacobian(),nor);
      }
      else
      {
         trans.SetIntPoint(&ip);
         
         w = ip.weight; // No detJ, cancels out for this integrator!
         w *= pcoeff->Eval(trans, ip); // p*w, evaluate thru face transformation
         if(coeff)
         {
            w *= coeff->Eval(trans,ip);
         }      

         // compute the normal (unnormalized!)
         CalcOrtho(trans.Jacobian(),nor);
      }     

      // shape nd x 1
      el.CalcShape(ip, shape);

      // f_f for each component
      for (int ii = 0; ii < dim; ii++)
      {
         temp_vec = shape; // nd x 1
         temp_vec *= w; // w*p * shape (nd x 1)
         temp_vec *= nor(ii); // w*p*n_ii * shape (nd x 1)
         elvect.AddSubVector(temp_vec,ii*nd);
      }
   }
}

void FollowerPressureLoad::AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &trans,
                                     const Vector &elfun,
                                     DenseMatrix &elmat)
{
   const int nd = el.GetDof();
   const int dim = el.GetDim() + 1;

   Vector shape, nor(dim);
   Vector Nk_dksi(dim);

   DenseMatrix Jacobian(dim), H1(dim), H2(dim), dshape(nd,dim), K_ik(dim);
   Jacobian = 0.0;
   H1 = 0.0;
   H2 = 0.0;
   K_ik = 0.0;
   dshape = 0.0;

   shape.SetSize(nd);
   shape = 0.0;
   nor = 0.0;

   elmat.SetSize(nd*dim);
   elmat = 0.0;

   if (formulation_type == form_type::TotalLagrangian)
   {
      EF.UseExternalData(elfun.GetData(), nd, dim);
   }

   double w = 0.0;
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // int order = 3 * el1.GetOrder() + 3; // correct order? Same as stress term
      // int order = trans.OrderGrad(&el1) + 2*el1.GetOrder(); // correct order?
      int order = 3*el.GetOrder() + trans.OrderW(); // correct order?
      ir = &IntRules.Get(trans.GetGeometryType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {     
      const IntegrationPoint &ip = ir->IntPoint(i);

      if (formulation_type == form_type::TotalLagrangian)
      {
         BuildCurrentTransformation(trans,ip,&EF);
         // curr_trans->SetIntPoint(&ip); // this is done inside above function
         
         w = ip.weight; // No detJ, cancels out for this integrator!
         w *= pcoeff->Eval(*curr_trans, ip); // p*w, evaluate thru face transformation
         if(coeff)
         {
            w *= coeff->Eval(*curr_trans,ip);
         }      

         // compute jacobian
         Jacobian = curr_trans->Jacobian();
      }
      else
      {
         trans.SetIntPoint(&ip);
         
         w = ip.weight; // No detJ, cancels out for this integrator!
         w *= pcoeff->Eval(trans, ip); // p*w, evaluate thru face transformation
         if(coeff)
         {
            w *= coeff->Eval(trans,ip);
         }

         // compute jacobian
         Jacobian = trans.Jacobian();   
      }

      // shape nd x 1
      el.CalcShape(ip, shape);

      // dshape 
      el.CalcDShape(ip,dshape); // nd x dim

      // setup H matrices
      // H1
      H1(0,1) = Jacobian(2,0);
      H1(1,0) = -H1(0,1);

      H1(0,2) = -Jacobian(1,0);
      H1(2,0) = -H1(0,2);

      H1(1,2) = Jacobian(0,0);
      H1(2,1) = -H1(1,2);

      // H2
      H2(0,1) = Jacobian(2,1);
      H2(1,0) = -H2(0,1);

      H2(0,2) = -Jacobian(1,1);
      H2(2,0) = -H2(0,2);

      H2(1,2) = Jacobian(0,1);
      H2(2,1) = -H2(1,2);
      
      for (int si = 0; si < nd; si++)
      {
         for (int sk = 0; sk < nd; sk++)
         {
            dshape.GetRow(sk,Nk_dksi);

            elmat(si,sk) += (Nk_dksi(0)*H2(0,0) - Nk_dksi(1)*H1(0,0)) * w * shape(si);
            elmat(si,sk + nd) += (Nk_dksi(0)*H2(0,1) - Nk_dksi(1)*H1(0,1)) * w * shape(si);
            if (dim == 3) elmat(si, sk + 2*nd) += (Nk_dksi(0)*H2(0,2) - Nk_dksi(1)*H1(0,2)) * w * shape(si);

            elmat(si + nd,sk) += (Nk_dksi(0)*H2(1,0) - Nk_dksi(1)*H1(1,0)) * w * shape(si);
            elmat(si + nd,sk + nd) += (Nk_dksi(0)*H2(1,1) - Nk_dksi(1)*H1(1,1)) * w * shape(si);
            if (dim == 3) elmat(si + nd, sk + 2*nd) += (Nk_dksi(0)*H2(1,2) - Nk_dksi(1)*H1(1,2)) * w * shape(si);

            elmat(si + 2*nd,sk) += (Nk_dksi(0)*H2(2,0) - Nk_dksi(1)*H1(2,0)) * w * shape(si);
            elmat(si + 2*nd,sk + nd) += (Nk_dksi(0)*H2(2,1) - Nk_dksi(1)*H1(2,1)) * w * shape(si);
            if (dim == 3) elmat(si + 2*nd, sk + 2*nd) += (Nk_dksi(0)*H2(2,2) - Nk_dksi(1)*H1(2,2)) * w * shape(si);
         }
      }
   }
}



// DG TERMS
// virtual external work for follower pressure load
void FollowerPressureLoad::AssembleFaceVector(const FiniteElement &el1, const FiniteElement &el2,
                                   FaceElementTransformations &Tr, const Vector &elfun, Vector &elvect)
{
   if (formulation_type == form_type::TotalLagrangian) MFEM_ABORT("TO DO");

   const int nd = el1.GetDof();
   const int dim = el1.GetDim();

   Vector shape, nor(dim), temp_vec(nd);

   shape.SetSize(nd);
   shape = 0.0;
   nor = 0.0;
   temp_vec = 0.0;

   elvect.SetSize(nd*dim);
   elvect = 0.0;

   double w = 0.0;
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order = 2 * el1.GetOrder(); // correct order? Same as vectorboundaryintegrator
      // int order = Tr.OrderGrad(&el1) + el1.GetOrder(); // correct order?
      // ir = &IntRules.Get(el1.GetGeomType(), order); // this is wrong, should be int rule for triangle not tet!!!!!
      ir = &IntRules.Get(Tr.GetGeometryType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {     
      const IntegrationPoint &ip = ir->IntPoint(i);
      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      // shape nd x 1
      el1.CalcShape(eip, shape);

      w = ip.weight; // No detJ, cancels out for this integrator!
      w *= pcoeff->Eval(Tr, ip); // p*w, evaluate thru face transformation
      if(coeff)
      {
         w *= coeff->Eval(Tr,ip);
      }      

      // compute the normal (unnormalized!)
      CalcOrtho(Tr.Jacobian(),nor);
      // // normalize nor, see example 18 in general it is not a unit normal
      // double normag = 0;
      // for (int ii = 0; ii < nor.Size(); ii++)
      // {
      //    normag += nor(ii)*nor(ii);
      // }
      // normag = sqrt(normag);
      // nor *= 1/normag;

      // f_f for each component
      for (int ii = 0; ii < dim; ii++)
      {
         temp_vec = shape; // nd x 1
         temp_vec *= w; // w*p * shape (nd x 1)
         temp_vec *= nor(ii); // w*p*n_ii * shape (nd x 1)
         elvect.AddSubVector(temp_vec,ii*nd);
      }
   }
}



// linearized virtual external work term for follower pressure load
void FollowerPressureLoad::AssembleFaceGrad(const FiniteElement &el1, const FiniteElement &el2,
                                   FaceElementTransformations &Tr, const Vector &elfun, DenseMatrix &elmat)
{
   const int nd = el1.GetDof();
   const int dim = el1.GetDim();

   Vector shape, nor(dim);

   DenseMatrix Jacobian(dim), H1(dim), H2(dim), dshape(nd,dim), K_ik(dim);
   Jacobian = 0.0;
   H1 = 0.0;
   H2 = 0.0;
   K_ik = 0.0;
   dshape = 0.0;

   shape.SetSize(nd);
   shape = 0.0;
   nor = 0.0;

   elmat.SetSize(nd*dim);
   elmat = 0.0;

   double w = 0.0;
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // int order = 3 * el1.GetOrder() + 3; // correct order? Same as stress term
      // int order = Tr.OrderGrad(&el1) + 2*el1.GetOrder(); // correct order?
      int order = 3*el1.GetOrder(); // correct order?
      // ir = &IntRules.Get(el1.GetGeomType(), order);
      ir = &IntRules.Get(Tr.GetGeometryType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {     
      const IntegrationPoint &ip = ir->IntPoint(i);
      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      // shape nd x 1
      el1.CalcShape(eip, shape);

      // dshape 
      el1.CalcDShape(eip,dshape); // nd x dim

      w = ip.weight; // No detJ, cancels out for this integrator!
      w *= pcoeff->Eval(Tr, ip); // p*w, evaluate thru face transformation
      if(coeff)
      {
         w *= coeff->Eval(Tr,ip);
      }
      if (curr_dt)
      {
         w *= alpha_f*beta*(*curr_dt)*(*curr_dt);
      }

      // compute jacobian
      Jacobian = Tr.Jacobian();

      // setup H matrices
      // H1
      H1(0,1) = Jacobian(2,0);
      H1(1,0) = -H1(0,1);

      H1(0,2) = -Jacobian(1,0);
      H1(2,0) = -H1(0,2);

      H1(1,2) = Jacobian(0,0);
      H1(2,1) = -H1(1,2);

      // H2
      H2(0,1) = Jacobian(2,1);
      H2(1,0) = -H2(0,1);

      H2(0,2) = -Jacobian(1,1);
      H2(2,0) = -H2(0,2);

      H2(1,2) = Jacobian(0,1);
      H2(2,1) = -H2(1,2);

      // form K_IK
      Vector Nk_dksi(dim);
      for (int si = 0; si < nd; si++)
      {
         for (int sk = 0; sk < nd; sk++)
         {
            K_ik = 0.0;

            dshape.GetRow(sk,Nk_dksi);

            // K_ik = NK_1*H2
            K_ik = H2;
            K_ik *= Nk_dksi(0);

            // K_ik = NK_1*H2 - NK_2*H1
            K_ik.AddMatrix(-Nk_dksi(1),H1,0,0);

            // K_ik = N_I * (NK_1*H2 - NK_2*H1)
            K_ik *= shape(si);

            // add into tangent stiffness matrix
            for (int d1 = 0; d1 < dim; d1++)
            {
               for (int d2 = 0; d2 < dim; d2++)
               {
                  elmat(si + nd*d1,sk + nd*d2) += w * K_ik(d1,d2); // here w = w * p
               }                  
            }
         }
      }
   }
}

void FollowerPressureLoad::BuildCurrentTransformation(ElementTransformation &T, const IntegrationPoint &ip, DenseMatrix *disp_pm)
{
   auto init_trans = dynamic_cast<IsoparametricTransformation *>(&T);
   auto elem = init_trans->GetFE();   

   if (curr_trans == nullptr || curr_trans->ElementNo != T.ElementNo)
   {
      if (curr_trans) delete curr_trans;
      curr_trans = new IsoparametricTransformation();
      curr_trans->SetFE(elem);
      
      DenseMatrix &init_pm = dynamic_cast<IsoparametricTransformation *>(&T)->GetPointMat();
      DenseMatrix curr_pm(init_pm);
      for (int ii = 0; ii < curr_pm.Height(); ii++)
      {
         curr_pm(0,ii) += (*disp_pm)(ii,0);
         curr_pm(1,ii) += (*disp_pm)(ii,1);
         curr_pm(2,ii) += (*disp_pm)(ii,2);
      }
      // curr_pm.Transpose(*disp_pm);

      // copy over T details
      curr_trans->SetPointMat(curr_pm);
      curr_trans->ElementNo = T.ElementNo;
      curr_trans->ElementType = T.ElementType;
   }
   // always set current int point
   curr_trans->SetIntPoint(&ip);
}