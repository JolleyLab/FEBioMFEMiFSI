#include "inc_constraint_form.hpp"

#include <iomanip>

using namespace mfem;
using namespace std;

// this integrator is always computed in TL approach
void IncConstraintIntegrator::AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &trans,
                                      const Vector &elfun,
                                      Vector &elvect)
{
   const int nd_p = el.GetDof();
   const int dim = el.GetDim();

   // number of nodes for displacement
   // if (fe_disp) delete fe_disp;
   auto fe_disp = u->FESpace()->GetFE(trans.ElementNo);
   const int nd_d = fe_disp->GetDof();

   Vector shape_p;
   DenseMatrix dshape_d, dshapedx, gradu, EF, F;

   shape_p.SetSize(nd_p);
   dshape_d.SetSize(nd_d, dim);
   dshapedx.SetSize(nd_d, dim);
   elvect.SetSize(nd_p);
   elvect = 0.0;
   gradu.SetSize(dim);
   gradu = 0.0;

   Vector vdof_vals(nd_d*dim);
   vdof_vals = 0.0;

   // get displacement dofs on the element
   u->GetElementDofValues(trans.ElementNo,vdof_vals); // double check this?
   F.SetSize(dim,dim);
   F = 0.0;
   EF.UseExternalData(vdof_vals.GetData(), nd_d, dim);

   double w;
   double Jminus1;
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order = el.GetOrder() + trans.OrderGrad(&el); // correct order?
      // int order = 2 * trans.OrderGrad(&el) + 3; // correct order?
      ir = &IntRules.Get(el.GetGeomType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {     
      const IntegrationPoint &ip = ir->IntPoint(i);
      trans.SetIntPoint(&ip);

      el.CalcShape(ip,shape_p);
      fe_disp->CalcDShape(ip,dshape_d);

      // ComputeDisplacementFE(ip, &el, dshape_d);

      // compute dshapedx using appropriate transformation for computing du/dX and F = dx/dX

      w = ip.weight; // get detJ for parent -> reference element always!
      if (formulation_type == form_type::TotalLagrangian)
      {
         Mult(dshape_d, trans.InverseJacobian(), dshapedx);
         w = trans.Weight();
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         BuildInitTransformation(trans,ip);
         Mult(dshape_d, init_trans->InverseJacobian(), dshapedx);
         w *= init_trans->Weight();
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }

      MultAtB(EF, dshapedx, gradu); // grad_phys u

      // compute detF - 1
      F = Identity_Mat;
      F += gradu;
      Jminus1 = F.Det();
      Jminus1 -= 1.0;

      for (int s = 0; s < nd_p; s++)
      {
         elvect(s) += w * shape_p(s) * Jminus1;
      }  
   }
   elvect *= -1.0;
}

void IncConstraintIntegrator::BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip)
{
   auto curr_trans = dynamic_cast<IsoparametricTransformation *>(&T);
   auto curr_elem = curr_trans->GetFE();

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
      Vector vals, x_orig(dim);
      x_orig = 0.0;
      Array<int> vdofs;

      // get vdofs depending on type of transformation
      if (T.ElementType == ElementTransformation::BDR_ELEMENT)
      {
         // boundary element transformation
         auto dof_trans = fespace->GetBdrElementVDofs(T.ElementNo,vdofs);
      }
      else if (T.ElementType == ElementTransformation::ELEMENT)
      {
         // element transformation
         auto dof_trans = fespace->GetElementVDofs(T.ElementNo,vdofs);
      }
      else
      {
         MFEM_ABORT("TO DO: ERROR, rebuilding init_transformation, not implemented element type");
      }
      
      orig_nodes->GetSubVector(vdofs,vals);
      // build point matrix for init_trans
      DenseMatrix init_pm(dim,curr_elem->GetDof());
      int idx, idx1, idx2;
      int vdim_num = vals.Size()/dim;
      int shift = vdofs.Size()/dim;
      for (int j=0; j<vdim_num; j++)
      {
         idx = j;
         idx1 = j + shift;
         if (dim == 3) { idx2 = j + shift*2; }

         x_orig(0) = vals(idx);
         x_orig(1) = vals(idx1);
         if (dim == 3) { x_orig(2) = vals(idx2); }

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

// void IncConstraintIntegrator::ComputeDisplacementFE(const IntegrationPoint &ip, FiniteElement *el, DenseMatrix &dshape, Vector &shape)
// {
//    // compute dshape and shape for displacement FE
//    int geom_type = el->GetGeomType();

//    if (fe == nullptr || geom_type != temp_geom_type)

//    if (geom_type == 2)
//    {
//       if(fe != nullptr) delete fe;
//       // triangle
//       fe = new H1_TriangleElement(orderd);
//    }
//    else if (geom_type == 3)
//    {
//       if(fe != nullptr) delete fe;
//       // square
//       fe = new H1_QuadrilateralElement(orderd);
//    }
//    else if (geom_type == 4)
//    {
//       if(fe != nullptr) delete fe;
//       // tet
//       fe = new H1_TetrahedronElement(orderd);
//    }
//    else if (geom_type == 5)
//    {
//       if(fe != nullptr) delete fe;
//       // cube
//       fe = new H1_HexahedronElement(orderd);
//    }
//    else if (geom_type == 6)
//    {
//       if(fe != nullptr) delete fe;
//       // prism
//       fe = new H1_WedgeElement(orderd);
//    }
//    else
//    {
//       MFEM_ABORT("Error: unsupported element type, see MFEM list");
//    }

//    // calc dshape and shape
//    fe->CalcShape(ip,shape);
//    fe->CalcDShape(ip,dshape);
// }