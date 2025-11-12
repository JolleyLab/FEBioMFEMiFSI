#include "ImmersedSolidIntegrator.hpp"

using namespace mfem;
using namespace std;
using namespace immersed_solid_serial;

void Immersed_Solid_Integrator::AssembleElementVector(const FiniteElement &el,
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

   inertial_force.SetSize(dim);
   bf_vec.SetSize(dim);

   // viscous fluid stress
   v_gf->FESpace()->GetElementVDofs(trans.ElementNo,v_vdofs);
   v_nodes.SetSize(nd*dim);
   v_gf->GetSubVector(v_vdofs,v_nodes);
   VF.UseExternalData(v_nodes.GetData(), nd, dim);
   dNd_currx.SetSize(nd,dim);
   e_v.SetSize(dim);
   Finv.SetSize(dim);
   dNidx_j.SetSize(dim);
   e_dNdx.SetSize(dim);

   if (formulation_type == form_type::TotalLagrangian)
   {
      F.SetSize(dim);
   }   

   // voigt
   Ni_dx.SetSize(dim);

   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
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

      // compute fluid viscous stress
      Finv = Identity_Mat;
      Finv += gradu;
      F = Finv;
      Finv.Invert();
      // dNd_currx
      Mult(dshapedx,Finv,dNd_currx);

      // epsilon * 2 * mu_f
      MultAtB(VF, dNd_currx, e_v); // dvdx
      e_v.Symmetrize();
      e_v *= 2.0*mu_f;
      if (subtract_fluid_visc_stress == false) e_v *= 0.0;

      // compute density
      double rhos = 0.0; // always variable used regardless of formulation
      if (formulation_type == form_type::TotalLagrangian)
      {
         rhos = rhos0_coeff->Eval(trans,ip);
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         double rhos0 = rhos0_coeff->Eval(*init_trans,ip); // eval on initial mesh
         rhos = rhos0 / mat_type->GetJ();
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }

      // subtract off fluid density
      // rho_diff = rhos_0 - rhof * J = J*(rhos - rhof)
      rhos -= formulation_type == form_type::TotalLagrangian ? rhof * mat_type->GetJ() : rhof;

      // acceleration
      a_gf->GetVectorValue(trans,ip,inertial_force);

      // body force (assumes forcing is applied over entire domain)
      for (auto &bf_term : *bf_terms)
      {
         if (formulation_type == form_type::TotalLagrangian)
         {
            bf_term.coeff->Eval(bf_vec,trans,ip); // eval on reference configuration!
         }
         else if (formulation_type == form_type::UpdatedLagrangian)
         {
            bf_term.coeff->Eval(bf_vec,*init_trans,ip); // eval on reference configuration!
         }
         else
         {
            MFEM_ABORT("Error: Unsupported formulation type");
         }         
         inertial_force.Add(-1.0,bf_vec);
      }

      // multiply by density
      // interial force = (rhos_0 - rhof J) * (d^2u/dt - f) = J*(rhos - rhof)*(d^2u/dt - f)
      inertial_force *= rhos;

      // compute element vector (voigt)
      if (formulation_type == form_type::UpdatedLagrangian)
      {
         for (int si = 0; si < nd; si++)
         {
            dNd_currx.GetRow(si,Ni_dx);
            // fluid terms
            e_v.Mult(Ni_dx,e_dNdx);
            
            Nix = Ni_dx(0);
            Niy = Ni_dx(1);
            if(dim == 3) Niz = Ni_dx(2);

            if (dim == 2)
            {
               elvect(si)      += w * (stress(0,0) * Nix + stress(0,1) * Niy) + w*shape(si)*inertial_force(0) - w * e_dNdx(0);
               elvect(si + nd) += w * (stress(1,1) * Niy + stress(0,1) * Nix) + w*shape(si)*inertial_force(1) - w * e_dNdx(1);
            }
            else if (dim == 3)
            {
               elvect(si)        += w * (stress(0,0) * Nix + stress(0,1) * Niy + stress(0,2) * Niz) + w*shape(si)*inertial_force(0) - w * e_dNdx(0);
               elvect(si + nd)   += w * (stress(1,1) * Niy + stress(0,1) * Nix + stress(1,2) * Niz) + w*shape(si)*inertial_force(1) - w * e_dNdx(1);
               elvect(si + 2*nd) += w * (stress(2,2) * Niz + stress(1,2) * Niy + stress(0,2) * Nix) + w*shape(si)*inertial_force(2) - w * e_dNdx(2);
            }            
         }
      }   
      else if(formulation_type == form_type::TotalLagrangian)
      {
         for (int si = 0; si < nd; si++)
         {
            // fluid terms
            dNd_currx.GetRow(si,dNidx_j);
            e_v.Mult(dNidx_j,e_dNdx);

            dshapedx.GetRow(si,Ni_dx);
            
            Nix = Ni_dx(0);
            Niy = Ni_dx(1);
            if(dim == 3) Niz = Ni_dx(2);

            if (dim == 2)
            {
               elvect(si)      += w * (stress(0,0) * F(0,0) * Nix + stress(1,1) * F(0,1) * Niy + stress(0,1)*( F(0,0)*Niy + F(0,1)*Nix ) ) + w*shape(si)*inertial_force(0) - w * mat_type->GetJ() * e_dNdx(0);
               elvect(si + nd) += w * (stress(0,0) * F(1,0) * Nix + stress(1,1) * F(1,1) * Niy + stress(0,1)*( F(1,0)*Niy + F(1,1)*Nix ) ) + w*shape(si)*inertial_force(1) - w * mat_type->GetJ() * e_dNdx(1);
            }
            else if (dim == 3)
            {
               elvect(si)        += w * (stress(0,0) * F(0,0) * Nix + stress(1,1) * F(0,1) * Niy + stress(2,2) * F(0,2) * Niz + stress(0,1)*( F(0,0)*Niy + F(0,1)*Nix ) + stress(1,2)*( F(0,1)*Niz + F(0,2)*Niy ) + stress(0,2)*( F(0,0)*Niz + F(0,2)*Nix ) ) + w*shape(si)*inertial_force(0) - w * mat_type->GetJ() * e_dNdx(0);
               elvect(si + nd)   += w * (stress(0,0) * F(1,0) * Nix + stress(1,1) * F(1,1) * Niy + stress(2,2) * F(1,2) * Niz + stress(0,1)*( F(1,0)*Niy + F(1,1)*Nix ) + stress(1,2)*( F(1,1)*Niz + F(1,2)*Niy ) + stress(0,2)*( F(1,0)*Niz + F(1,2)*Nix ) ) + w*shape(si)*inertial_force(1) - w * mat_type->GetJ() * e_dNdx(1);
               elvect(si + 2*nd) += w * (stress(0,0) * F(2,0) * Nix + stress(1,1) * F(2,1) * Niy + stress(2,2) * F(2,2) * Niz + stress(0,1)*( F(2,0)*Niy + F(2,1)*Nix ) + stress(1,2)*( F(2,1)*Niz + F(2,2)*Niy ) + stress(0,2)*( F(2,0)*Niz + F(2,2)*Nix ) ) + w*shape(si)*inertial_force(2) - w * mat_type->GetJ() * e_dNdx(2);
            }            
         }
      }
   }

   // mfem::out << "Element id = " << trans.ElementNo << std::endl;
   // auto Nodes = el.GetNodes();
   // mfem::Vector curr_pos(3);
   // for(int j=0; j<Nodes.GetNPoints(); ++j)
   // {
   //    auto ip_node = Nodes.IntPoint(j);
   //    trans.SetIntPoint(&ip_node);
   //    trans.Transform(ip_node, curr_pos);
   //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
   // }
   // mfem::out << "Element Vector = " << std::endl;
   // for (int ii=0; ii<nd; ++ii)
   // {
   //    mfem::out << elvect(ii) << " ";
   //    mfem::out << elvect(ii + nd) << " ";
   //    if (dim == 3) mfem::out << elvect(ii + 2*nd) << " ";
   //    mfem::out << std::endl;
   // }
   // // mfem::out << "VDOFS = " << std::endl;
   // // v_vdofs.Print();
   // // for (int ii=0; ii<nd; ++ii)
   // // {
   // //    mfem::out << "Node " << ii << " = (" << v_vdofs[ii] << ", " << v_vdofs[ii + nd] << ", " << v_vdofs[ii + 2*nd] << ")" << std::endl;
   // // }
   // // elvect.Print();
   // double temp = 0.0;
}



// generic linearized virtual internal work term
void Immersed_Solid_Integrator::AssembleElementGrad(const FiniteElement &el,
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

   // additional terms 
   Ni_dx.SetSize(dim);
   Nk_dx.SetSize(dim);
   temp_vec.SetSize(dim);
   gradphys_shape.SetSize(nd,dim);

   // fluid additional terms
   dNdx_i.SetSize(dim);
   dNdx_k.SetSize(dim);
   FinvT_dNdX.SetSize(dim);

   // viscous fluid stress
   v_gf->FESpace()->GetElementVDofs(trans.ElementNo,v_vdofs);
   v_nodes.SetSize(nd*dim);
   v_gf->GetSubVector(v_vdofs,v_nodes);
   VF.UseExternalData(v_nodes.GetData(), nd, dim);
   dNd_currx.SetSize(nd,dim);
   e_v.SetSize(dim);
   Finv.SetSize(dim);
   dNidx_j.SetSize(dim);
   e_dNdx.SetSize(dim);

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
      int order = 2 * el.GetOrder() + trans.OrderW();
      ir = &IntRules.Get(el.GetGeomType(), order);
   }

   for (int i = 0; i < ir->GetNPoints(); i++)
   {      
      const IntegrationPoint &ip = ir->IntPoint(i);
      trans.SetIntPoint(&ip);

      el.CalcShape(ip, shape);
      el.CalcDShape(ip, dshape);

      // compute dshapedx (dNdX) using appropriate transformation for computing du/dX and F = dx/dX
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
      
      // compute fluid viscous stress
      Finv = Identity_Mat;
      Finv += gradu;
      F = Finv;
      Finv.Invert();
      // dNd_currx
      Mult(dshapedx,Finv,dNd_currx); // dNdx = dN/dX * dX/dx = dN/dX * F^{-1}

      // epsilon * 2 * mu_f
      MultAtB(VF, dNd_currx, e_v); // dvdx
      e_v.Symmetrize();
      e_v *= 2.0*mu_f;
      if (subtract_fluid_visc_stress == false) e_v *= 0.0;

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

      // compute density
      double rhos = 0.0; // always variable used regardless of formulation
      if (formulation_type == form_type::TotalLagrangian)
      {
         rhos = rhos0_coeff->Eval(trans,ip);
      }
      else if (formulation_type == form_type::UpdatedLagrangian)
      {
         double rhos0 = rhos0_coeff->Eval(*init_trans,ip); // eval on initial mesh
         rhos = rhos0 / mat_type->GetJ();
      }
      else
      {
         MFEM_ABORT("Error: Unsupported formulation type");
      }

      // subtract off fluid density
      // rho_diff = rhos_0 - rhof * J = J*(rhos - rhof)
      rhos -= formulation_type == form_type::TotalLagrangian ? rhof * mat_type->GetJ() : rhof;

      // acceleration
      a_gf->GetVectorValue(trans,ip,inertial_force);

      // body force (assumes forcing is applied over entire domain)
      for (auto &bf_term : *bf_terms)
      {
         if (formulation_type == form_type::TotalLagrangian)
         {
            bf_term.coeff->Eval(bf_vec,trans,ip); // eval on reference configuration!
         }
         else if (formulation_type == form_type::UpdatedLagrangian)
         {
            bf_term.coeff->Eval(bf_vec,*init_trans,ip); // eval on reference configuration!
         }
         else
         {
            MFEM_ABORT("Error: Unsupported formulation type");
         }     
         inertial_force.Add(-1.0,bf_vec);
      }

      double am = alpha_m;
      double av = alpha_f*gamma*(*dt);
      double af = alpha_f*beta*(*dt)*(*dt);

      double inertial_val = 0.0;
      if (formulation_type == form_type::UpdatedLagrangian)
      {
         double fac;
         double fac2;

         double int_val;
         double visc_val;

         FinvT = F;
         FinvT.Invert();
         FinvT.Transpose();

         // element stiffness //
         for (int si = 0; si < nd; si++)
         {
            dshapedx.GetRow(si,Ni_dx); // this is dNdX
            dNd_currx.GetRow(si,dNdx_i); // this is dNdx
            e_v.Mult(dNdx_i,e_dNdx);
            
            Nix = dNdx_i(0);
            Niy = dNdx_i(1);
            if(dim == 3) Niz = dNdx_i(2);

            for (int sk = 0; sk < nd; sk++)
            {
               dshapedx.GetRow(sk,Nk_dx); // this is dNdX
               FinvT.Mult(Nk_dx,FinvT_dNdX); 
               dNd_currx.GetRow(sk,dNdx_k); // this is dNdx

               Nkx = dNdx_k(0);
               Nky = dNdx_k(1);
               if(dim == 3) Nkz = dNdx_k(2);

               ////// dR/dA //////
               // INERTIAL STIFFNESS //
               inertial_val = w*am*rhos*shape(si)*shape(sk);
               elmat(si,sk) += inertial_val;
               elmat(si + nd, sk + nd) += inertial_val;
               if (dim == 3) elmat(si + 2*nd, sk + 2*nd) += inertial_val;

               ////// dR/dV //////
               // VELOCITY TERM FROM FLUID //
               if (subtract_fluid_visc_stress)
               {
                  fac = -w*mu_f* av;
                  fac2 = dNdx_i*dNdx_k;
                  if (dim == 2)
                  {
                     elmat(si,sk)            += fac * (fac2 + dNd_currx(si,0)*dNd_currx(sk,0));
                     elmat(si,sk + nd)       += fac * dNd_currx(si,1)*dNd_currx(sk,0);
                     elmat(si + nd,sk)       += fac * dNd_currx(si,0)*dNd_currx(sk,1);
                     elmat(si + nd, sk + nd) += fac * (fac2 + dNd_currx(si,1)*dNd_currx(sk,1));
                  }
                  else if (dim == 3)
                  {
                     elmat(si,sk)             += fac * (fac2 + dNd_currx(si,0)*dNd_currx(sk,0));
                     elmat(si,sk+nd)          += fac * dNd_currx(si,1)*dNd_currx(sk,0);
                     elmat(si,sk+2*nd)        += fac * dNd_currx(si,2)*dNd_currx(sk,0);

                     elmat(si+nd,sk)          += fac * dNd_currx(si,0)*dNd_currx(sk,1);
                     elmat(si+nd,sk+nd)       += fac * (fac2 + dNd_currx(si,1)*dNd_currx(sk,1));
                     elmat(si+nd,sk+2*nd)     += fac * dNd_currx(si,2)*dNd_currx(sk,1);

                     elmat(si+2*nd,sk)        += fac * dNd_currx(si,0)*dNd_currx(sk,2);
                     elmat(si+2*nd,sk+nd)     += fac * dNd_currx(si,1)*dNd_currx(sk,2);
                     elmat(si+2*nd,sk+2*nd)   += fac * (fac2 + dNd_currx(si,2)*dNd_currx(sk,2));
                  }
               }

               ////// dR/dU //////
               // SHAPE DERIVATIVE TERMS //
               int_val = -af*w*shape(si)*rhof;
               visc_val = subtract_fluid_visc_stress ? -af*w : 0.0;
               if (dim == 2)
               {
                  elmat(si,sk)            += int_val*inertial_force(0)*FinvT_dNdX(0) + visc_val *e_dNdx(0)*FinvT_dNdX(0);
                  elmat(si,sk + nd)       += int_val*inertial_force(0)*FinvT_dNdX(1) + visc_val *e_dNdx(0)*FinvT_dNdX(1);
                  elmat(si + nd,sk)       += int_val*inertial_force(1)*FinvT_dNdX(0) + visc_val *e_dNdx(1)*FinvT_dNdX(0);
                  elmat(si + nd, sk + nd) += int_val*inertial_force(1)*FinvT_dNdX(1) + visc_val *e_dNdx(1)*FinvT_dNdX(1);
               }
               else if (dim == 3)
               {
                  elmat(si,sk)             += int_val*inertial_force(0)*FinvT_dNdX(0) + visc_val *e_dNdx(0)*FinvT_dNdX(0);
                  elmat(si,sk+nd)          += int_val*inertial_force(0)*FinvT_dNdX(1) + visc_val *e_dNdx(0)*FinvT_dNdX(1);
                  elmat(si,sk+2*nd)        += int_val*inertial_force(0)*FinvT_dNdX(2) + visc_val *e_dNdx(0)*FinvT_dNdX(2);

                  elmat(si+nd,sk)          += int_val*inertial_force(1)*FinvT_dNdX(0) + visc_val *e_dNdx(1)*FinvT_dNdX(0);
                  elmat(si+nd,sk+nd)       += int_val*inertial_force(1)*FinvT_dNdX(1) + visc_val *e_dNdx(1)*FinvT_dNdX(1);
                  elmat(si+nd,sk+2*nd)     += int_val*inertial_force(1)*FinvT_dNdX(2) + visc_val *e_dNdx(1)*FinvT_dNdX(2);

                  elmat(si+2*nd,sk)        += int_val*inertial_force(2)*FinvT_dNdX(0) + visc_val *e_dNdx(2)*FinvT_dNdX(0);
                  elmat(si+2*nd,sk+nd)     += int_val*inertial_force(2)*FinvT_dNdX(1) + visc_val *e_dNdx(2)*FinvT_dNdX(1);
                  elmat(si+2*nd,sk+2*nd)   += int_val*inertial_force(2)*FinvT_dNdX(2) + visc_val *e_dNdx(2)*FinvT_dNdX(2);
               }

               // MATERIAL STIFFNESS //
               AssembleMaterialStiffness(si,sk,dim,nd,af,elmat);  

               // GEOMETRIC STIFFNESS //
               stress.Mult(dNdx_k,temp_vec);
               geo_val = dNdx_i * temp_vec * w * af;               

               elmat(si,sk) += geo_val;
               elmat(si + nd, sk + nd) += geo_val;
               if (dim == 3) elmat(si + 2*nd, sk + 2*nd) += geo_val;
            }
         }
      }
      else if (formulation_type == form_type::TotalLagrangian)
      {
         // element stiffness //
         FinvT = F;
         FinvT.Invert();
         FinvT.Transpose();

         double fac;
         double fac2;

         double int_val;
         double visc_val;

         for (int si = 0; si < nd; si++)
         {
            dshapedx.GetRow(si,Ni_dx); // this is dNdX
            dNd_currx.GetRow(si,dNdx_i); // this is dNdx
            e_v.Mult(dNdx_i,e_dNdx);
            
            Nix = Ni_dx(0);
            Niy = Ni_dx(1);
            if(dim == 3) Niz = Ni_dx(2);

            for (int sk = 0; sk < nd; sk++)
            {
               dshapedx.GetRow(sk,Nk_dx); // this is dNdX
               FinvT.Mult(Nk_dx,FinvT_dNdX); // this is dNdx
               dNd_currx.GetRow(sk,dNdx_k);

               Nkx = Nk_dx(0);
               Nky = Nk_dx(1);
               if(dim == 3) Nkz = Nk_dx(2);

               ////// dR/dA //////
               // INERTIAL STIFFNESS //
               inertial_val = w*am*rhos*shape(si)*shape(sk);
               elmat(si,sk) += inertial_val;
               elmat(si + nd, sk + nd) += inertial_val;
               if (dim == 3) elmat(si + 2*nd, sk + 2*nd) += inertial_val;

               ////// dR/dV //////
               // VELOCITY TERM FROM FLUID //
               if (subtract_fluid_visc_stress)
               {
                  fac = -w*mu_f*mat_type->GetJ() * av;
                  fac2 = dNdx_i*dNdx_k;
                  if (dim == 2)
                  {
                     elmat(si,sk)            += fac * (fac2 + dNd_currx(si,0)*dNd_currx(sk,0));
                     elmat(si,sk + nd)       += fac * dNd_currx(si,1)*dNd_currx(sk,0);
                     elmat(si + nd,sk)       += fac * dNd_currx(si,0)*dNd_currx(sk,1);
                     elmat(si + nd, sk + nd) += fac * (fac2 + dNd_currx(si,1)*dNd_currx(sk,1));
                  }
                  else if (dim == 3)
                  {
                     elmat(si,sk)             += fac * (fac2 + dNd_currx(si,0)*dNd_currx(sk,0));
                     elmat(si,sk+nd)          += fac * dNd_currx(si,1)*dNd_currx(sk,0);
                     elmat(si,sk+2*nd)        += fac * dNd_currx(si,2)*dNd_currx(sk,0);

                     elmat(si+nd,sk)          += fac * dNd_currx(si,0)*dNd_currx(sk,1);
                     elmat(si+nd,sk+nd)       += fac * (fac2 + dNd_currx(si,1)*dNd_currx(sk,1));
                     elmat(si+nd,sk+2*nd)     += fac * dNd_currx(si,2)*dNd_currx(sk,1);

                     elmat(si+2*nd,sk)        += fac * dNd_currx(si,0)*dNd_currx(sk,2);
                     elmat(si+2*nd,sk+nd)     += fac * dNd_currx(si,1)*dNd_currx(sk,2);
                     elmat(si+2*nd,sk+2*nd)   += fac * (fac2 + dNd_currx(si,2)*dNd_currx(sk,2));
                  }
               }

               ////// dR/dU //////
               // SHAPE DERIVATIVE TERMS //
               int_val = -af*w*mat_type->GetJ()*shape(si)*rhof;
               visc_val = subtract_fluid_visc_stress ? -af*w*mat_type->GetJ() : 0.0;
               if (dim == 2)
               {
                  elmat(si,sk)            += int_val*inertial_force(0)*FinvT_dNdX(0) + visc_val *e_dNdx(0)*FinvT_dNdX(0);
                  elmat(si,sk + nd)       += int_val*inertial_force(0)*FinvT_dNdX(1) + visc_val *e_dNdx(0)*FinvT_dNdX(1);
                  elmat(si + nd,sk)       += int_val*inertial_force(1)*FinvT_dNdX(0) + visc_val *e_dNdx(1)*FinvT_dNdX(0);
                  elmat(si + nd, sk + nd) += int_val*inertial_force(1)*FinvT_dNdX(1) + visc_val *e_dNdx(1)*FinvT_dNdX(1);
               }
               else if (dim == 3)
               {
                  elmat(si,sk)             += int_val*inertial_force(0)*FinvT_dNdX(0) + visc_val *e_dNdx(0)*FinvT_dNdX(0);
                  elmat(si,sk+nd)          += int_val*inertial_force(0)*FinvT_dNdX(1) + visc_val *e_dNdx(0)*FinvT_dNdX(1);
                  elmat(si,sk+2*nd)        += int_val*inertial_force(0)*FinvT_dNdX(2) + visc_val *e_dNdx(0)*FinvT_dNdX(2);

                  elmat(si+nd,sk)          += int_val*inertial_force(1)*FinvT_dNdX(0) + visc_val *e_dNdx(1)*FinvT_dNdX(0);
                  elmat(si+nd,sk+nd)       += int_val*inertial_force(1)*FinvT_dNdX(1) + visc_val *e_dNdx(1)*FinvT_dNdX(1);
                  elmat(si+nd,sk+2*nd)     += int_val*inertial_force(1)*FinvT_dNdX(2) + visc_val *e_dNdx(1)*FinvT_dNdX(2);

                  elmat(si+2*nd,sk)        += int_val*inertial_force(2)*FinvT_dNdX(0) + visc_val *e_dNdx(2)*FinvT_dNdX(0);
                  elmat(si+2*nd,sk+nd)     += int_val*inertial_force(2)*FinvT_dNdX(1) + visc_val *e_dNdx(2)*FinvT_dNdX(1);
                  elmat(si+2*nd,sk+2*nd)   += int_val*inertial_force(2)*FinvT_dNdX(2) + visc_val *e_dNdx(2)*FinvT_dNdX(2);
               }

               // INTERNAL WORK TERMS//
               // MATERIAL STIFFNESS //
               AssembleMaterialStiffness(si,sk,dim,nd,af,elmat); 

               // GEOMETRIC STIFFNESS //
               stress.Mult(Nk_dx,temp_vec);
               geo_val = Ni_dx * temp_vec * w * af;

               elmat(si,sk) += geo_val;
               elmat(si + nd, sk + nd) += geo_val;
               if (dim == 3) elmat(si + 2*nd, sk + 2*nd) += geo_val;
            }
         }
      }      
   }

   // mfem::out << "\n";
   // mfem::out << "Element id = " << trans.ElementNo << std::endl;
   // auto Nodes = el.GetNodes();
   // mfem::Vector curr_pos(3);
   // for(int j=0; j<Nodes.GetNPoints(); ++j)
   // {
   //    auto ip_node = Nodes.IntPoint(j);
   //    trans.SetIntPoint(&ip_node);
   //    trans.Transform(ip_node, curr_pos);
   //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
   // }

   // // Open a file to write the matrix
   // std::ofstream outfile("matrix_ke_mfem.txt", std::ios::out);

   // // Print the matrix ke to the file
   // Array<int> febio_map(3*nd);
   // for (int i=0; i<nd; ++i)
   // {
   //    febio_map[i] = 3*i;
   //    febio_map[i+nd] = 3*i+1;
   //    if (dim == 3) febio_map[i+2*nd] = 3*i+2;
   // }

   // // int ii,jj;
   // // for (int i = 0; i < 3*nd; ++i) 
   // // {
   // //    ii = febio_map[i];
   // //    for (int j = 0; j < 3*nd; ++j) 
   // //    {
   // //       jj = febio_map[j];
   // //       if (abs(elmat(ii,jj)) >= 1e-11) outfile << std::setw(10) << std::setprecision(4) << elmat(ii,jj);
   // //       else outfile << std::setw(10) << std::setprecision(4) << 0.0;
   // //       if (j < 3*nd - 1) outfile << " ";
   // //    }
   // //    outfile << std::endl;
   // // }
   // // outfile << std::endl; // Add an extra line between matrices

   // int ii,jj;
   // int precision_out = 6;
   // int width_out = precision_out + 6;
   // for (int i = 0; i < 3*nd; ++i) 
   // {
   //    ii = i;
   //    for (int j = 0; j < 3*nd; ++j) 
   //    {
   //       jj = j;
   //       if (abs(elmat(ii,jj)) >= 1e-12) outfile << std::setw(width_out) << std::setprecision(precision_out) << elmat(ii,jj);
   //       else outfile << std::setw(width_out) << std::setprecision(precision_out) << 0.0;
   //       if (j < 3*nd - 1) outfile << " ";
   //    }
   //    outfile << std::endl;
   // }
   // outfile << std::endl; // Add an extra line between matrices


   // // Close the file
   // outfile.close();
   // int temp_val = 0;
}

void Immersed_Solid_Integrator::AssembleMaterialStiffness(int &si,int &sk, const int &dim, const int &nd, double &af, DenseMatrix &elmat)
{
   if (formulation_type == form_type::UpdatedLagrangian)
   {
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
         elmat(si,sk)            += (Nix * CBK(0,0) + Niy * CBK(2,0)) * w * af;
         elmat(si,sk + nd)       += (Nix * CBK(0,1) + Niy * CBK(2,1)) * w * af;
         elmat(si + nd,sk)       += (Niy * CBK(1,0) + Nix * CBK(2,0)) * w * af;
         elmat(si + nd, sk + nd) += (Niy * CBK(1,1) + Nix * CBK(2,1)) * w * af;
      }
      else if (dim == 3)
      {
         elmat(si,sk)             += (Nix*CBK(0,0) + Niy*CBK(3,0) + Niz*CBK(5,0)) * w * af;
         elmat(si,sk+nd)          += (Nix*CBK(0,1) + Niy*CBK(3,1) + Niz*CBK(5,1)) * w * af;
         elmat(si,sk+2*nd)        += (Nix*CBK(0,2) + Niy*CBK(3,2) + Niz*CBK(5,2)) * w * af;

         elmat(si+nd,sk)          += (Niy*CBK(1,0) + Nix*CBK(3,0) + Niz*CBK(4,0)) * w * af;
         elmat(si+nd,sk+nd)       += (Niy*CBK(1,1) + Nix*CBK(3,1) + Niz*CBK(4,1)) * w * af;
         elmat(si+nd,sk+2*nd)     += (Niy*CBK(1,2) + Nix*CBK(3,2) + Niz*CBK(4,2)) * w * af;

         elmat(si+2*nd,sk)        += (Niz*CBK(2,0) + Niy*CBK(4,0) + Nix*CBK(5,0)) * w * af;
         elmat(si+2*nd,sk+nd)     += (Niz*CBK(2,1) + Niy*CBK(4,1) + Nix*CBK(5,1)) * w * af;
         elmat(si+2*nd,sk+2*nd)   += (Niz*CBK(2,2) + Niy*CBK(4,2) + Nix*CBK(5,2)) * w * af;
      }
   }
   else if (formulation_type == form_type::TotalLagrangian)
   {
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
         elmat(si,sk)            += (F(0,0) * Nix * CBK(0,0) + F(0,1) * Niy * CBK(1,0) + ( F(0,0)*Niy + F(0,1)*Nix ) * CBK(2,0)) * w * af;
         elmat(si,sk + nd)       += (F(0,0) * Nix * CBK(0,1) + F(0,1) * Niy * CBK(1,1) + ( F(0,0)*Niy + F(0,1)*Nix ) * CBK(2,1)) * w * af;

         elmat(si + nd,sk)       += (F(1,0) * Nix * CBK(0,0) + F(1,1) * Niy * CBK(1,0) + ( F(1,0)*Niy + F(1,1)*Nix ) * CBK(2,0)) * w * af;
         elmat(si + nd, sk + nd) += (F(1,0) * Nix * CBK(0,1) + F(1,1) * Niy * CBK(1,1) + ( F(1,0)*Niy + F(1,1)*Nix ) * CBK(2,1)) * w * af;               
      }
      else if (dim == 3)
      {
         elmat(si,sk)             += (CBK(0,0) * F(0,0) * Nix + CBK(1,0) * F(0,1) * Niy + CBK(2,0) * F(0,2) * Niz + CBK(3,0)*( F(0,0)*Niy + F(0,1)*Nix ) + CBK(4,0)*( F(0,1)*Niz + F(0,2)*Niy ) + CBK(5,0)*( F(0,0)*Niz + F(0,2)*Nix )) * w * af;
         elmat(si,sk+nd)          += (CBK(0,1) * F(0,0) * Nix + CBK(1,1) * F(0,1) * Niy + CBK(2,1) * F(0,2) * Niz + CBK(3,1)*( F(0,0)*Niy + F(0,1)*Nix ) + CBK(4,1)*( F(0,1)*Niz + F(0,2)*Niy ) + CBK(5,1)*( F(0,0)*Niz + F(0,2)*Nix )) * w * af;
         elmat(si,sk+2*nd)        += (CBK(0,2) * F(0,0) * Nix + CBK(1,2) * F(0,1) * Niy + CBK(2,2) * F(0,2) * Niz + CBK(3,2)*( F(0,0)*Niy + F(0,1)*Nix ) + CBK(4,2)*( F(0,1)*Niz + F(0,2)*Niy ) + CBK(5,2)*( F(0,0)*Niz + F(0,2)*Nix )) * w * af;

         elmat(si+nd,sk)          += (CBK(0,0) * F(1,0) * Nix + CBK(1,0) * F(1,1) * Niy + CBK(2,0) * F(1,2) * Niz + CBK(3,0)*( F(1,0)*Niy + F(1,1)*Nix ) + CBK(4,0)*( F(1,1)*Niz + F(1,2)*Niy ) + CBK(5,0)*( F(1,0)*Niz + F(1,2)*Nix )) * w * af;
         elmat(si+nd,sk+nd)       += (CBK(0,1) * F(1,0) * Nix + CBK(1,1) * F(1,1) * Niy + CBK(2,1) * F(1,2) * Niz + CBK(3,1)*( F(1,0)*Niy + F(1,1)*Nix ) + CBK(4,1)*( F(1,1)*Niz + F(1,2)*Niy ) + CBK(5,1)*( F(1,0)*Niz + F(1,2)*Nix )) * w * af;
         elmat(si+nd,sk+2*nd)     += (CBK(0,2) * F(1,0) * Nix + CBK(1,2) * F(1,1) * Niy + CBK(2,2) * F(1,2) * Niz + CBK(3,2)*( F(1,0)*Niy + F(1,1)*Nix ) + CBK(4,2)*( F(1,1)*Niz + F(1,2)*Niy ) + CBK(5,2)*( F(1,0)*Niz + F(1,2)*Nix )) * w * af;

         elmat(si+2*nd,sk)        += (CBK(0,0) * F(2,0) * Nix + CBK(1,0) * F(2,1) * Niy + CBK(2,0) * F(2,2) * Niz + CBK(3,0)*( F(2,0)*Niy + F(2,1)*Nix ) + CBK(4,0)*( F(2,1)*Niz + F(2,2)*Niy ) + CBK(5,0)*( F(2,0)*Niz + F(2,2)*Nix )) * w * af;
         elmat(si+2*nd,sk+nd)     += (CBK(0,1) * F(2,0) * Nix + CBK(1,1) * F(2,1) * Niy + CBK(2,1) * F(2,2) * Niz + CBK(3,1)*( F(2,0)*Niy + F(2,1)*Nix ) + CBK(4,1)*( F(2,1)*Niz + F(2,2)*Niy ) + CBK(5,1)*( F(2,0)*Niz + F(2,2)*Nix )) * w * af;
         elmat(si+2*nd,sk+2*nd)   += (CBK(0,2) * F(2,0) * Nix + CBK(1,2) * F(2,1) * Niy + CBK(2,2) * F(2,2) * Niz + CBK(3,2)*( F(2,0)*Niy + F(2,1)*Nix ) + CBK(4,2)*( F(2,1)*Niz + F(2,2)*Niy ) + CBK(5,2)*( F(2,0)*Niz + F(2,2)*Nix )) * w * af;
      }
   }     
}

void Immersed_Solid_Integrator::BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip)
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