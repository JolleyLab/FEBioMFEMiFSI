#include "StressTerm.hpp"

#include <iomanip>

using namespace mfem;
using namespace std;
using namespace mat;

// generic virtual internal work term
void StressIntegrator::AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &trans,
                                      const Vector &elfun,
                                      Vector &elvect)
{
   const int nd = el.GetDof();
   const int dim = el.GetDim();

   u.SetSize(dim);
   vec.SetSize(nd);

   shape.SetSize(nd);
   dshape.SetSize(nd, dim);
   dshapedx.SetSize(nd, dim);
   elvect.SetSize(nd * dim);
   elvect = 0.0;
   gradu.SetSize(dim);
   gradu = 0.0;
   stress.SetSize(dim);
   stress = 0.0;
   gradphys_shape.SetSize(nd,dim);

   EF.UseExternalData(elfun.GetData(), nd, dim);

   if (formulation_type == form_type::TotalLagrangian)
   {
      F.SetSize(dim);
   }   

   // voigt
   Ni_dx.SetSize(dim);

   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order = 2 * el.GetOrder() + trans.OrderW(); // correct order?
      // int order = 2 * trans.OrderGrad(&el) + trans.OrderW(); // correct order?
      // int order = 2 * trans.OrderGrad(&el) + 3; // correct order?
      // int order = 2 * el.GetOrder() + 3; // correct order?
      ir = &IntRules.Get(el.GetGeomType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {      
      const IntegrationPoint &ip = ir->IntPoint(i);
      trans.SetIntPoint(&ip);

      el.CalcShape(ip, shape);
      el.CalcDShape(ip, dshape);

      // compute dshapedx using appropriate transformation for computing du/dX and F = dx/dX
      if (formulation_type == form_type::TotalLagrangian)
      {
         Mult(dshape, trans.InverseJacobian(), dshapedx); // (dN/dX)
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         BuildInitTransformation(trans,ip);
         Mult(dshape, init_trans->InverseJacobian(), dshapedx); // (dN/dX)
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }

      w = ip.weight * trans.Weight();

      if (coeff)
      {
         w *= coeff->Eval(trans, ip);
      }

      MultAtB(EF, dshapedx, gradu); // dudX
      EF.MultTranspose(shape, u); // u
      
      // compute strain measures, material properties at the current integration point, stored for all operations later on
      if (formulation_type == form_type::TotalLagrangian)
      {
         mat_type->SetStrainMeasures(u,gradu,ip,trans);
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         mat_type->SetStrainMeasures(u,gradu,ip,*init_trans);
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }      

      if (incompressible_flag)
      {
         double p = pres->GetValue(*init_trans,ip);
         mat_type->SetPressure(p);
      }

      // compute stress
      if (formulation_type == form_type::TotalLagrangian)
      {
         mat_type->EvalPK2Stress(stress);
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         mat_type->EvalCauchyStress(stress);      
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }

      // compute element vector (voigt)
      if (formulation_type == form_type::UpdatedLagrangian)
      {
         Mult(dshape,trans.InverseJacobian(),gradphys_shape); // gradref_phi * J^{-1} nd x dim (dN/dx)

         for (int si = 0; si < nd; si++)
         {
            gradphys_shape.GetRow(si,Ni_dx);
            
            Nix = Ni_dx(0);
            Niy = Ni_dx(1);
            if(dim == 3) Niz = Ni_dx(2);

            if (dim == 2)
            {
               elvect(si)      += w * (stress(0,0) * Nix + stress(0,1) * Niy);
               elvect(si + nd) += w * (stress(1,1) * Niy + stress(0,1) * Nix);
            }
            else if (dim == 3)
            {
               elvect(si)        += w * (stress(0,0) * Nix + stress(0,1) * Niy + stress(0,2) * Niz);
               elvect(si + nd)   += w * (stress(1,1) * Niy + stress(0,1) * Nix + stress(1,2) * Niz);
               elvect(si + 2*nd) += w * (stress(2,2) * Niz + stress(1,2) * Niy + stress(0,2) * Nix);
            }            
         }
      }   
      else if(formulation_type == form_type::TotalLagrangian)
      {
         F = Identity_Mat;
         F += gradu;

         for (int si = 0; si < nd; si++)
         {
            dshapedx.GetRow(si,Ni_dx);
            
            Nix = Ni_dx(0);
            Niy = Ni_dx(1);
            if(dim == 3) Niz = Ni_dx(2);

            if (dim == 2)
            {
               elvect(si)      += w * (stress(0,0) * F(0,0) * Nix + stress(1,1) * F(0,1) * Niy + stress(0,1)*( F(0,0)*Niy + F(0,1)*Nix ) );
               elvect(si + nd) += w * (stress(0,0) * F(1,0) * Nix + stress(1,1) * F(1,1) * Niy + stress(0,1)*( F(1,0)*Niy + F(1,1)*Nix ) );
            }
            else if (dim == 3)
            {
               elvect(si)        += w * (stress(0,0) * F(0,0) * Nix + stress(1,1) * F(0,1) * Niy + stress(2,2) * F(0,2) * Niz + stress(0,1)*( F(0,0)*Niy + F(0,1)*Nix ) + stress(1,2)*( F(0,1)*Niz + F(0,2)*Niy ) + stress(0,2)*( F(0,0)*Niz + F(0,2)*Nix ) );
               elvect(si + nd)   += w * (stress(0,0) * F(1,0) * Nix + stress(1,1) * F(1,1) * Niy + stress(2,2) * F(1,2) * Niz + stress(0,1)*( F(1,0)*Niy + F(1,1)*Nix ) + stress(1,2)*( F(1,1)*Niz + F(1,2)*Niy ) + stress(0,2)*( F(1,0)*Niz + F(1,2)*Nix ) );
               elvect(si + 2*nd) += w * (stress(0,0) * F(2,0) * Nix + stress(1,1) * F(2,1) * Niy + stress(2,2) * F(2,2) * Niz + stress(0,1)*( F(2,0)*Niy + F(2,1)*Nix ) + stress(1,2)*( F(2,1)*Niz + F(2,2)*Niy ) + stress(0,2)*( F(2,0)*Niz + F(2,2)*Nix ) );
            }            
         }
      }
      
      // std::cout << "------------\n";
      // std::cout << setprecision(17);
      // elvect.Print(cout);    
   }
   // std::cout << setprecision(17);
   // std::cout << "elem no =" << trans.ElementNo << std::endl;
   // elvect.Print(cout);
}



// generic linearized virtual internal work term
void StressIntegrator::AssembleElementGrad(const FiniteElement &el,
                                     ElementTransformation &trans,
                                     const Vector &elfun,
                                     DenseMatrix &elmat)
{
   const int nd = el.GetDof();
   const int dim = el.GetDim();

   shape.SetSize(nd);
   shape = 0.0;
   u.SetSize(dim);
   u = 0.0;
   dshape.SetSize(nd, dim);
   dshape = 0.0;
   dshapedx.SetSize(nd, dim);
   dshapedx = 0.0;
   gradu.SetSize(dim);
   gradu = 0.0;
   stress.SetSize(dim);
   stress = 0.0;
   
   EF.UseExternalData(elfun.GetData(), nd, dim);

   elmat.SetSize(nd*dim, nd*dim);
   elmat = 0.0;

   // DenseMatrix geo_mat(dim), mat_mat(dim);
   // geo_mat = 0.0;
   // mat_mat = 0.0;

   // additional terms 
   Ni_dx.SetSize(dim);
   Nk_dx.SetSize(dim);
   temp_vec.SetSize(dim);
   gradphys_shape.SetSize(nd,dim);

   if (formulation_type == form_type::TotalLagrangian)
   {
      F.SetSize(dim);
   }

   // voigt terms
   if (dim == 2)
   {
      CBK.SetSize(3,2);
   }
   else if (dim == 3)
   {
      CBK.SetSize(6,3);
   }

   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // // int order = 3 * trans.OrderGrad(&el) + trans.OrderW(); // correct order?
      // int order = 3 * el.GetOrder() + trans.OrderW(); // correct order?
      // // int order = 2 * el.GetOrder() + 3; // correct order?

      // old
      // int order = 3 * el.GetOrder() + trans.OrderW(); // correct order?

      // 5-1-24
      int order = 2 * el.GetOrder() + trans.OrderW();
      ir = &IntRules.Get(el.GetGeomType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {      
      const IntegrationPoint &ip = ir->IntPoint(i);
      trans.SetIntPoint(&ip);

      el.CalcShape(ip, shape);
      el.CalcDShape(ip, dshape);

      // compute dshapedx using appropriate transformation for computing du/dX and F = dx/dX
      if (formulation_type == form_type::TotalLagrangian)
      {
         Mult(dshape, trans.InverseJacobian(), dshapedx); // gradrefN * J^{-1} // nd x dim
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         BuildInitTransformation(trans,ip);
         Mult(dshape, init_trans->InverseJacobian(), dshapedx); // gradrefN * J^{-1} // nd x dim
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }
      
      w = ip.weight * trans.Weight(); // w * detJ for integration

      if (coeff)
      {
         w *= coeff->Eval(trans, ip);
      }

      MultAtB(EF, dshapedx, gradu); // grad_phys u
      EF.MultTranspose(shape, u); // u

      // compute strain measures, material properties at the current integration point, stored for all operations later on
      if (formulation_type == form_type::TotalLagrangian)
      {
         mat_type->SetStrainMeasures(u,gradu,ip,trans);
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         mat_type->SetStrainMeasures(u,gradu,ip,*init_trans);
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }

      if (incompressible_flag)
      {
         double p = pres->GetValue(*init_trans,ip);
         mat_type->SetPressure(p);
      }

      // compute stress
      if (formulation_type == form_type::TotalLagrangian)
      {
         mat_type->EvalPK2Stress(stress);
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         mat_type->EvalCauchyStress(stress);
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }

      // get voigt elasticity tensor
      if (dim == 3)
      {
         if (formulation_type == form_type::UpdatedLagrangian)
         {
            D = mat_type->GetElasticityTensorVoigt();
         }
         else
         {
            D = mat_type->GetElasticityTensorVoigt();
         }         
      }

      if (formulation_type == form_type::UpdatedLagrangian)
      {
         // N,x
         Mult(dshape,trans.InverseJacobian(),gradphys_shape); // gradref_phi * J^{-1} // nd x dim

         // element stiffness //
         for (int si = 0; si < nd; si++)
         {
            gradphys_shape.GetRow(si,Ni_dx); // this is dNdx
            
            Nix = Ni_dx(0);
            Niy = Ni_dx(1);
            if(dim == 3) Niz = Ni_dx(2);

            for (int sk = 0; sk < nd; sk++)
            {
               gradphys_shape.GetRow(sk,Nk_dx); 

               Nkx = Nk_dx(0);
               Nky = Nk_dx(1);
               if(dim == 3) Nkz = Nk_dx(2);

               // MATERIAL STIFFNESS //
               // [C]BK
               if (dim == 2)
               {
                  CBK(0,0) = mat_type->EvalElasticityTensor(0,0,0,0) * Nkx + mat_type->EvalElasticityTensor(0,0,0,1) * Nky;
                  CBK(0,1) = mat_type->EvalElasticityTensor(0,0,1,1) * Nky + mat_type->EvalElasticityTensor(0,0,0,1) * Nkx;

                  CBK(1,0) = mat_type->EvalElasticityTensor(1,1,0,0) * Nkx + mat_type->EvalElasticityTensor(1,1,0,1) * Nky;
                  CBK(1,1) = mat_type->EvalElasticityTensor(1,1,1,1) * Nky + mat_type->EvalElasticityTensor(1,1,0,1) * Nkx;

                  CBK(2,0) = mat_type->EvalElasticityTensor(0,1,0,0) * Nkx + mat_type->EvalElasticityTensor(0,1,0,1) * Nky;
                  CBK(2,1) = mat_type->EvalElasticityTensor(0,1,1,1) * Nky + mat_type->EvalElasticityTensor(0,1,0,1) * Nkx;
               }
               else if (dim == 3)
               {
                  CBK(0,0) = D->Elem(0,0) * Nkx + D->Elem(0,3) * Nky + D->Elem(0,5) * Nkz;
                  CBK(0,1) = D->Elem(0,1) * Nky + D->Elem(0,3) * Nkx + D->Elem(0,4) * Nkz;
                  CBK(0,2) = D->Elem(0,2) * Nkz + D->Elem(0,4) * Nky + D->Elem(0,5) * Nkx;

                  CBK(1,0) = D->Elem(1,0) * Nkx + D->Elem(1,3) * Nky + D->Elem(1,5) * Nkz;
                  CBK(1,1) = D->Elem(1,1) * Nky + D->Elem(1,3) * Nkx + D->Elem(1,4) * Nkz;
                  CBK(1,2) = D->Elem(1,2) * Nkz + D->Elem(1,4) * Nky + D->Elem(1,5) * Nkx;

                  CBK(2,0) = D->Elem(2,0) * Nkx + D->Elem(2,3) * Nky + D->Elem(2,5) * Nkz;
                  CBK(2,1) = D->Elem(2,1) * Nky + D->Elem(2,3) * Nkx + D->Elem(2,4) * Nkz;
                  CBK(2,2) = D->Elem(2,2) * Nkz + D->Elem(2,4) * Nky + D->Elem(2,5) * Nkx;

                  CBK(3,0) = D->Elem(3,0) * Nkx + D->Elem(3,3) * Nky + D->Elem(3,5) * Nkz;
                  CBK(3,1) = D->Elem(3,1) * Nky + D->Elem(3,3) * Nkx + D->Elem(3,4) * Nkz;
                  CBK(3,2) = D->Elem(3,2) * Nkz + D->Elem(3,4) * Nky + D->Elem(3,5) * Nkx;

                  CBK(4,0) = D->Elem(4,0) * Nkx + D->Elem(4,3) * Nky + D->Elem(4,5) * Nkz;
                  CBK(4,1) = D->Elem(4,1) * Nky + D->Elem(4,3) * Nkx + D->Elem(4,4) * Nkz;
                  CBK(4,2) = D->Elem(4,2) * Nkz + D->Elem(4,4) * Nky + D->Elem(4,5) * Nkx;

                  CBK(5,0) = D->Elem(5,0) * Nkx + D->Elem(5,3) * Nky + D->Elem(5,5) * Nkz;
                  CBK(5,1) = D->Elem(5,1) * Nky + D->Elem(5,3) * Nkx + D->Elem(5,4) * Nkz;
                  CBK(5,2) = D->Elem(5,2) * Nkz + D->Elem(5,4) * Nky + D->Elem(5,5) * Nkx;
               }     

               // BI^T [C] BK
               if (dim == 2)
               {
                  elmat(si,sk)            += (Nix * CBK(0,0) + Niy * CBK(2,0)) * w;
                  elmat(si,sk + nd)       += (Nix * CBK(0,1) + Niy * CBK(2,1)) * w;
                  elmat(si + nd,sk)       += (Niy * CBK(1,0) + Nix * CBK(2,0)) * w;
                  elmat(si + nd, sk + nd) += (Niy * CBK(1,1) + Nix * CBK(2,1)) * w;
               }
               else if (dim == 3)
               {
                  elmat(si,sk)             += (Nix*CBK(0,0) + Niy*CBK(3,0) + Niz*CBK(5,0)) * w;
                  elmat(si,sk+nd)          += (Nix*CBK(0,1) + Niy*CBK(3,1) + Niz*CBK(5,1)) * w;
                  elmat(si,sk+2*nd)        += (Nix*CBK(0,2) + Niy*CBK(3,2) + Niz*CBK(5,2)) * w;

                  elmat(si+nd,sk)          += (Niy*CBK(1,0) + Nix*CBK(3,0) + Niz*CBK(4,0)) * w;
                  elmat(si+nd,sk+nd)       += (Niy*CBK(1,1) + Nix*CBK(3,1) + Niz*CBK(4,1)) * w;
                  elmat(si+nd,sk+2*nd)     += (Niy*CBK(1,2) + Nix*CBK(3,2) + Niz*CBK(4,2)) * w;

                  elmat(si+2*nd,sk)        += (Niz*CBK(2,0) + Niy*CBK(4,0) + Nix*CBK(5,0)) * w;
                  elmat(si+2*nd,sk+nd)     += (Niz*CBK(2,1) + Niy*CBK(4,1) + Nix*CBK(5,1)) * w;
                  elmat(si+2*nd,sk+2*nd)   += (Niz*CBK(2,2) + Niy*CBK(4,2) + Nix*CBK(5,2)) * w;
               }     

               // GEOMETRIC STIFFNESS //
               stress.Mult(Nk_dx,temp_vec);
               geo_val = Ni_dx * temp_vec * w;

               elmat(si,sk) += geo_val;
               elmat(si + nd, sk + nd) += geo_val;
               if (dim == 3) elmat(si + 2*nd, sk + 2*nd) += geo_val;
            }
         }
      }
      else if (formulation_type == form_type::TotalLagrangian)
      {
         // element stiffness //
         F = Identity_Mat;
         F += gradu;

         for (int si = 0; si < nd; si++)
         {
            dshapedx.GetRow(si,Ni_dx); // this is dNdX
            
            Nix = Ni_dx(0);
            Niy = Ni_dx(1);
            if(dim == 3) Niz = Ni_dx(2);

            for (int sk = 0; sk < nd; sk++)
            {
               dshapedx.GetRow(sk,Nk_dx);

               Nkx = Nk_dx(0);
               Nky = Nk_dx(1);
               if(dim == 3) Nkz = Nk_dx(2);

               // MATERIAL STIFFNESS //
               // [C]BK
               if (dim == 2)
               {
                  CBK(0,0) = mat_type->EvalElasticityTensor(0,0,0,0) * F(0,0) * Nkx + mat_type->EvalElasticityTensor(0,0,1,1) * F(0,1) * Nky + mat_type->EvalElasticityTensor(0,0,0,1)*( F(0,0)*Nky + F(0,1)*Nkx );
                  CBK(0,1) = mat_type->EvalElasticityTensor(0,0,0,0) * F(1,0) * Nkx + mat_type->EvalElasticityTensor(0,0,1,1) * F(1,1) * Nky + mat_type->EvalElasticityTensor(0,0,0,1)*( F(1,0)*Nky + F(1,1)*Nkx );

                  CBK(1,0) = mat_type->EvalElasticityTensor(1,1,0,0) * F(0,0) * Nkx + mat_type->EvalElasticityTensor(1,1,1,1) * F(0,1) * Nky + mat_type->EvalElasticityTensor(1,1,0,1)*( F(0,0)*Nky + F(0,1)*Nkx );
                  CBK(1,1) = mat_type->EvalElasticityTensor(1,1,0,0) * F(1,0) * Nkx + mat_type->EvalElasticityTensor(1,1,1,1) * F(1,1) * Nky + mat_type->EvalElasticityTensor(1,1,0,1)*( F(1,0)*Nky + F(1,1)*Nkx );

                  CBK(2,0) = mat_type->EvalElasticityTensor(0,1,0,0) * F(0,0) * Nkx + mat_type->EvalElasticityTensor(0,1,1,1) * F(0,1) * Nky + mat_type->EvalElasticityTensor(0,1,0,1)*( F(0,0)*Nky + F(0,1)*Nkx );
                  CBK(2,1) = mat_type->EvalElasticityTensor(0,1,0,0) * F(1,0) * Nkx + mat_type->EvalElasticityTensor(0,1,1,1) * F(1,1) * Nky + mat_type->EvalElasticityTensor(0,1,0,1)*( F(1,0)*Nky + F(1,1)*Nkx );
               }
               else if (dim == 3)
               {
                  CBK(0,0) = D->Elem(0,0) * F(0,0) * Nkx + D->Elem(0,1) * F(0,1) * Nky + D->Elem(0,2) * F(0,2) * Nkz + D->Elem(0,3)*( F(0,0)*Nky + F(0,1)*Nkx ) + D->Elem(0,4)*( F(0,1)*Nkz + F(0,2)*Nky ) + D->Elem(0,5)*( F(0,0)*Nkz + F(0,2)*Nkx );
                  CBK(0,1) = D->Elem(0,0) * F(1,0) * Nkx + D->Elem(0,1) * F(1,1) * Nky + D->Elem(0,2) * F(1,2) * Nkz + D->Elem(0,3)*( F(1,0)*Nky + F(1,1)*Nkx ) + D->Elem(0,4)*( F(1,1)*Nkz + F(1,2)*Nky ) + D->Elem(0,5)*( F(1,0)*Nkz + F(1,2)*Nkx );
                  CBK(0,2) = D->Elem(0,0) * F(2,0) * Nkx + D->Elem(0,1) * F(2,1) * Nky + D->Elem(0,2) * F(2,2) * Nkz + D->Elem(0,3)*( F(2,0)*Nky + F(2,1)*Nkx ) + D->Elem(0,4)*( F(2,1)*Nkz + F(2,2)*Nky ) + D->Elem(0,5)*( F(2,0)*Nkz + F(2,2)*Nkx );

                  CBK(1,0) = D->Elem(1,0) * F(0,0) * Nkx + D->Elem(1,1) * F(0,1) * Nky + D->Elem(1,2) * F(0,2) * Nkz + D->Elem(1,3)*( F(0,0)*Nky + F(0,1)*Nkx ) + D->Elem(1,4)*( F(0,1)*Nkz + F(0,2)*Nky ) + D->Elem(1,5)*( F(0,0)*Nkz + F(0,2)*Nkx );
                  CBK(1,1) = D->Elem(1,0) * F(1,0) * Nkx + D->Elem(1,1) * F(1,1) * Nky + D->Elem(1,2) * F(1,2) * Nkz + D->Elem(1,3)*( F(1,0)*Nky + F(1,1)*Nkx ) + D->Elem(1,4)*( F(1,1)*Nkz + F(1,2)*Nky ) + D->Elem(1,5)*( F(1,0)*Nkz + F(1,2)*Nkx );
                  CBK(1,2) = D->Elem(1,0) * F(2,0) * Nkx + D->Elem(1,1) * F(2,1) * Nky + D->Elem(1,2) * F(2,2) * Nkz + D->Elem(1,3)*( F(2,0)*Nky + F(2,1)*Nkx ) + D->Elem(1,4)*( F(2,1)*Nkz + F(2,2)*Nky ) + D->Elem(1,5)*( F(2,0)*Nkz + F(2,2)*Nkx );

                  CBK(2,0) = D->Elem(2,0) * F(0,0) * Nkx + D->Elem(2,1) * F(0,1) * Nky + D->Elem(2,2) * F(0,2) * Nkz + D->Elem(2,3)*( F(0,0)*Nky + F(0,1)*Nkx ) + D->Elem(2,4)*( F(0,1)*Nkz + F(0,2)*Nky ) + D->Elem(2,5)*( F(0,0)*Nkz + F(0,2)*Nkx );
                  CBK(2,1) = D->Elem(2,0) * F(1,0) * Nkx + D->Elem(2,1) * F(1,1) * Nky + D->Elem(2,2) * F(1,2) * Nkz + D->Elem(2,3)*( F(1,0)*Nky + F(1,1)*Nkx ) + D->Elem(2,4)*( F(1,1)*Nkz + F(1,2)*Nky ) + D->Elem(2,5)*( F(1,0)*Nkz + F(1,2)*Nkx );
                  CBK(2,2) = D->Elem(2,0) * F(2,0) * Nkx + D->Elem(2,1) * F(2,1) * Nky + D->Elem(2,2) * F(2,2) * Nkz + D->Elem(2,3)*( F(2,0)*Nky + F(2,1)*Nkx ) + D->Elem(2,4)*( F(2,1)*Nkz + F(2,2)*Nky ) + D->Elem(2,5)*( F(2,0)*Nkz + F(2,2)*Nkx );

                  CBK(3,0) = D->Elem(3,0) * F(0,0) * Nkx + D->Elem(3,1) * F(0,1) * Nky + D->Elem(3,2) * F(0,2) * Nkz + D->Elem(3,3)*( F(0,0)*Nky + F(0,1)*Nkx ) + D->Elem(3,4)*( F(0,1)*Nkz + F(0,2)*Nky ) + D->Elem(3,5)*( F(0,0)*Nkz + F(0,2)*Nkx );
                  CBK(3,1) = D->Elem(3,0) * F(1,0) * Nkx + D->Elem(3,1) * F(1,1) * Nky + D->Elem(3,2) * F(1,2) * Nkz + D->Elem(3,3)*( F(1,0)*Nky + F(1,1)*Nkx ) + D->Elem(3,4)*( F(1,1)*Nkz + F(1,2)*Nky ) + D->Elem(3,5)*( F(1,0)*Nkz + F(1,2)*Nkx );
                  CBK(3,2) = D->Elem(3,0) * F(2,0) * Nkx + D->Elem(3,1) * F(2,1) * Nky + D->Elem(3,2) * F(2,2) * Nkz + D->Elem(3,3)*( F(2,0)*Nky + F(2,1)*Nkx ) + D->Elem(3,4)*( F(2,1)*Nkz + F(2,2)*Nky ) + D->Elem(3,5)*( F(2,0)*Nkz + F(2,2)*Nkx );

                  CBK(4,0) = D->Elem(4,0) * F(0,0) * Nkx + D->Elem(4,1) * F(0,1) * Nky + D->Elem(4,2) * F(0,2) * Nkz + D->Elem(4,3)*( F(0,0)*Nky + F(0,1)*Nkx ) + D->Elem(4,4)*( F(0,1)*Nkz + F(0,2)*Nky ) + D->Elem(4,5)*( F(0,0)*Nkz + F(0,2)*Nkx );
                  CBK(4,1) = D->Elem(4,0) * F(1,0) * Nkx + D->Elem(4,1) * F(1,1) * Nky + D->Elem(4,2) * F(1,2) * Nkz + D->Elem(4,3)*( F(1,0)*Nky + F(1,1)*Nkx ) + D->Elem(4,4)*( F(1,1)*Nkz + F(1,2)*Nky ) + D->Elem(4,5)*( F(1,0)*Nkz + F(1,2)*Nkx );
                  CBK(4,2) = D->Elem(4,0) * F(2,0) * Nkx + D->Elem(4,1) * F(2,1) * Nky + D->Elem(4,2) * F(2,2) * Nkz + D->Elem(4,3)*( F(2,0)*Nky + F(2,1)*Nkx ) + D->Elem(4,4)*( F(2,1)*Nkz + F(2,2)*Nky ) + D->Elem(4,5)*( F(2,0)*Nkz + F(2,2)*Nkx );

                  CBK(5,0) = D->Elem(5,0) * F(0,0) * Nkx + D->Elem(5,1) * F(0,1) * Nky + D->Elem(5,2) * F(0,2) * Nkz + D->Elem(5,3)*( F(0,0)*Nky + F(0,1)*Nkx ) + D->Elem(5,4)*( F(0,1)*Nkz + F(0,2)*Nky ) + D->Elem(5,5)*( F(0,0)*Nkz + F(0,2)*Nkx );
                  CBK(5,1) = D->Elem(5,0) * F(1,0) * Nkx + D->Elem(5,1) * F(1,1) * Nky + D->Elem(5,2) * F(1,2) * Nkz + D->Elem(5,3)*( F(1,0)*Nky + F(1,1)*Nkx ) + D->Elem(5,4)*( F(1,1)*Nkz + F(1,2)*Nky ) + D->Elem(5,5)*( F(1,0)*Nkz + F(1,2)*Nkx );
                  CBK(5,2) = D->Elem(5,0) * F(2,0) * Nkx + D->Elem(5,1) * F(2,1) * Nky + D->Elem(5,2) * F(2,2) * Nkz + D->Elem(5,3)*( F(2,0)*Nky + F(2,1)*Nkx ) + D->Elem(5,4)*( F(2,1)*Nkz + F(2,2)*Nky ) + D->Elem(5,5)*( F(2,0)*Nkz + F(2,2)*Nkx );
               }     

               // BI^T [C] BK
               if (dim == 2)
               {
                  elmat(si,sk)            += (F(0,0) * Nix * CBK(0,0) + F(0,1) * Niy * CBK(1,0) + ( F(0,0)*Niy + F(0,1)*Nix ) * CBK(2,0)) * w;
                  elmat(si,sk + nd)       += (F(0,0) * Nix * CBK(0,1) + F(0,1) * Niy * CBK(1,1) + ( F(0,0)*Niy + F(0,1)*Nix ) * CBK(2,1)) * w;

                  elmat(si + nd,sk)       += (F(1,0) * Nix * CBK(0,0) + F(1,1) * Niy * CBK(1,0) + ( F(1,0)*Niy + F(1,1)*Nix ) * CBK(2,0)) * w;
                  elmat(si + nd, sk + nd) += (F(1,0) * Nix * CBK(0,1) + F(1,1) * Niy * CBK(1,1) + ( F(1,0)*Niy + F(1,1)*Nix ) * CBK(2,1)) * w;
               }
               else if (dim == 3)
               {
                  elmat(si,sk)             += (CBK(0,0) * F(0,0) * Nix + CBK(1,0) * F(0,1) * Niy + CBK(2,0) * F(0,2) * Niz + CBK(3,0)*( F(0,0)*Niy + F(0,1)*Nix ) + CBK(4,0)*( F(0,1)*Niz + F(0,2)*Niy ) + CBK(5,0)*( F(0,0)*Niz + F(0,2)*Nix )) * w;
                  elmat(si,sk+nd)          += (CBK(0,1) * F(0,0) * Nix + CBK(1,1) * F(0,1) * Niy + CBK(2,1) * F(0,2) * Niz + CBK(3,1)*( F(0,0)*Niy + F(0,1)*Nix ) + CBK(4,1)*( F(0,1)*Niz + F(0,2)*Niy ) + CBK(5,1)*( F(0,0)*Niz + F(0,2)*Nix )) * w;
                  elmat(si,sk+2*nd)        += (CBK(0,2) * F(0,0) * Nix + CBK(1,2) * F(0,1) * Niy + CBK(2,2) * F(0,2) * Niz + CBK(3,2)*( F(0,0)*Niy + F(0,1)*Nix ) + CBK(4,2)*( F(0,1)*Niz + F(0,2)*Niy ) + CBK(5,2)*( F(0,0)*Niz + F(0,2)*Nix )) * w;

                  elmat(si+nd,sk)          += (CBK(0,0) * F(1,0) * Nix + CBK(1,0) * F(1,1) * Niy + CBK(2,0) * F(1,2) * Niz + CBK(3,0)*( F(1,0)*Niy + F(1,1)*Nix ) + CBK(4,0)*( F(1,1)*Niz + F(1,2)*Niy ) + CBK(5,0)*( F(1,0)*Niz + F(1,2)*Nix )) * w;
                  elmat(si+nd,sk+nd)       += (CBK(0,1) * F(1,0) * Nix + CBK(1,1) * F(1,1) * Niy + CBK(2,1) * F(1,2) * Niz + CBK(3,1)*( F(1,0)*Niy + F(1,1)*Nix ) + CBK(4,1)*( F(1,1)*Niz + F(1,2)*Niy ) + CBK(5,1)*( F(1,0)*Niz + F(1,2)*Nix )) * w;
                  elmat(si+nd,sk+2*nd)     += (CBK(0,2) * F(1,0) * Nix + CBK(1,2) * F(1,1) * Niy + CBK(2,2) * F(1,2) * Niz + CBK(3,2)*( F(1,0)*Niy + F(1,1)*Nix ) + CBK(4,2)*( F(1,1)*Niz + F(1,2)*Niy ) + CBK(5,2)*( F(1,0)*Niz + F(1,2)*Nix )) * w;

                  elmat(si+2*nd,sk)        += (CBK(0,0) * F(2,0) * Nix + CBK(1,0) * F(2,1) * Niy + CBK(2,0) * F(2,2) * Niz + CBK(3,0)*( F(2,0)*Niy + F(2,1)*Nix ) + CBK(4,0)*( F(2,1)*Niz + F(2,2)*Niy ) + CBK(5,0)*( F(2,0)*Niz + F(2,2)*Nix )) * w;
                  elmat(si+2*nd,sk+nd)     += (CBK(0,1) * F(2,0) * Nix + CBK(1,1) * F(2,1) * Niy + CBK(2,1) * F(2,2) * Niz + CBK(3,1)*( F(2,0)*Niy + F(2,1)*Nix ) + CBK(4,1)*( F(2,1)*Niz + F(2,2)*Niy ) + CBK(5,1)*( F(2,0)*Niz + F(2,2)*Nix )) * w;
                  elmat(si+2*nd,sk+2*nd)   += (CBK(0,2) * F(2,0) * Nix + CBK(1,2) * F(2,1) * Niy + CBK(2,2) * F(2,2) * Niz + CBK(3,2)*( F(2,0)*Niy + F(2,1)*Nix ) + CBK(4,2)*( F(2,1)*Niz + F(2,2)*Niy ) + CBK(5,2)*( F(2,0)*Niz + F(2,2)*Nix )) * w;
               }     

               // GEOMETRIC STIFFNESS //
               stress.Mult(Nk_dx,temp_vec);
               geo_val = Ni_dx * temp_vec * w;

               elmat(si,sk) += geo_val;
               elmat(si + nd, sk + nd) += geo_val;
               if (dim == 3) elmat(si + 2*nd, sk + 2*nd) += geo_val;
            }
         }
      }      
   }
   // std::cout << "------------------------\n";
   // std::cout << "Elem num = " << trans.ElementNo << std::endl;
   // std::cout << "geo mat = " << std::endl;
   // geo_mat.PrintMatlab(cout);
   // std::cout << "------------------------\n";
   // std::cout << "mat mat = " << std::endl;
   // mat_mat.PrintMatlab(cout);
   // std::cout << "------------------------\n";
}

void StressIntegrator::BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip)
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
      int nd = vals.Size()/dim;
      int shift = vdofs.Size()/dim;
      for (int j=0; j<nd; j++)
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