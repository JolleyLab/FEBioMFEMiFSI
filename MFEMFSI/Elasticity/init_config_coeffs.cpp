#include "init_config_coeffs.hpp"

void InitConfigDispCoeff::Eval(Vector &u, ElementTransformation &T, const IntegrationPoint &ip)
{      
   // build initial transformation
   BuildInitTransformation(T,ip);
   
   // evaluate total disp ON INITIAL MESH
   Vector td(dim);
   td = 0.0;
   total_disp->Eval(td,*init_trans,ip);
   
   // set displacement
   u.SetSize(dim);
   u = 0.0;
   u = td;
};

void InitConfigDispCoeff::BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip)
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

void InitConfigDispCoeff::SetTime(double t_)
{
   time = t_;
   total_disp->SetTime(time);
}

// scalar version
double InitConfigScalarCoeff::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{      
   // build initial transformation
   BuildInitTransformation(T,ip);

   // evaluate total disp on INITIAL MESH
   return total_disp_comp->Eval(*init_trans,ip);
};

void InitConfigScalarCoeff::BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip)
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

void InitConfigScalarCoeff::SetTime(double t_)
{
   time = t_;
   total_disp_comp->SetTime(time);
}