#include "IncompressibleFluidIntegrator.hpp"
#include <general/forall.hpp>

#include "../../Utilities/linear_element_ext.hpp"

using namespace mfem;
using namespace std;
using namespace ifs_solver;

/////////////////////////////////////////////////////////////////////////////////////

void ifs_solver::IncompressibleFluidIntegrator::AssembleElementVector(const Array<const FiniteElement *> &el,
                                       ElementTransformation &trans,
                                       const Array<const Vector *> &elfun,
                                       const Array<Vector *> &elvec)
{
   if (el.Size() != 2)
   {
      mfem_error("Finite element space of incorrect block number");
   }

   ///////// setup data /////////
   int nd_u = el[0]->GetDof();
   int nd_p = el[1]->GetDof();

   int dim = el[0]->GetDim();
   int spaceDim = trans.GetSpaceDim();

   if (dim != spaceDim)
   {
      mfem_error("Not defined on manifold meshes");
   }

   dshape_u.SetSize(nd_u, dim);
   dshapedx_u.SetSize(nd_u, dim);
   shape_u.SetSize(nd_u);
   gradu.SetSize(dim);
   u.SetSize(dim);
   u_dot_gradu.SetSize(dim);
   EF.UseExternalData(elfun[0]->GetData(),nd_u,dim); // velocity dofs
   shape_u_lap.SetSize(nd_u);
   shape_u_lap = 0.0;
   Nadx_sigma.SetSize(dim);
   u_lap.SetSize(dim);
   u_lap = 0.0;
   if (mesh_velo) u_mesh.SetSize(dim);

   dshape_p.SetSize(nd_p, dim);
   shape_p.SetSize(nd_p);
   dshapedx_p.SetSize(nd_p, dim);
   gradp.SetSize(dim);

   rm.SetSize(dim);
   rm_inertial.SetSize(dim);
   Nadx.SetSize(dim);
   gradu_rm.SetSize(dim);

   fluid_stress.SetSize(dim);
   f_vec.SetSize(dim);
   
   double w = 0.0;
   elvec[0]->SetSize(nd_u*dim);
   (*elvec[0]) = 0.0;
   elvec[1]->SetSize(nd_p);
   (*elvec[1]) = 0.0;

   ///////// integration rule /////////
   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      int order = 2*el[0]->GetOrder(); // exact integration for linear terms
      ir = &IntRules.Get(el[0]->GetGeomType(), order);
   }

   ///////// loop over integration points /////////
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      ///////// compute data /////////
      const IntegrationPoint &ip = ir->IntPoint(i);
      trans.SetIntPoint(&ip);

      // pressure shape functions
      el[1]->CalcShape(ip,shape_p);
      p = (*elfun[1]) * shape_p;
      el[1]->CalcDShape(ip, dshape_p);
      Mult(dshape_p, trans.InverseJacobian(), dshapedx_p);
      dshapedx_p.MultTranspose(*elfun[1],gradp); // grad p

      // velocity shape functions
      el[0]->CalcShape(ip,shape_u);
      el[0]->CalcDShape(ip,dshape_u);
      Mult(dshape_u, trans.InverseJacobian(), dshapedx_u);
      MultAtB(EF, dshapedx_u, gradu); // grad u
      EF.MultTranspose(shape_u,u); // u
      if (mesh_velo)
      {
         mesh_velo->Eval(u_mesh,trans,ip);
         u -= u_mesh; // u = (u - u_mesh)
      }
      gradu.Mult(u,u_dot_gradu); // u dot gradu
      divu = gradu.Trace();   
      
      if (nxx_flag)
      {
         // Not needed for p=1, minimal effect if included
         // this will be zero on linear triangles and tets
         // H1FECollection: Will work for all element types except wedges and pyramids
         el[0]->CalcPhysLaplacian(trans,shape_u_lap);
         EF.MultTranspose(shape_u_lap,u_lap); // u_i,jj
      } 

      ///// Compute rm at the integration point /////
      // rm = rho * (udot + u*gradu - f) + gradp - mu*d2u/dxdx
      // rm_inertial = rho * (udot + u*gradu - f)
      double mu_f = dyn_vis->Eval(trans,ip);
      Compute_rm(rm,rm_inertial,mu_f,trans,ip,el);

      // fluid stress
      fluid_stress.Set(2.0*mu_f,gradu);
      fluid_stress.Symmetrize();
      fluid_stress(0,0) -= p;
      fluid_stress(1,1) -= p;
      if (dim == 3) fluid_stress(2,2) -= p;      

      // stabilization parameters
      if (supg_flag) tau_m = tau_m_evaluator->EvalTauM(trans,ip,u);
      if (vms_flag) tau_b = tau_m_evaluator->EvalTauB(trans,ip,tau_m,rm);
      if (lsic_flag) tau_c = mass_stab_coeff->Eval_Coeff(tau_m);

      ///// Add contributions to element vector /////
      w = ip.weight * trans.Weight(); // w = w_q * detJ
      double rho = rhof->Eval(trans,ip);      

      ///////// momentum residual /////////
      double visc_valx = 0.0;
      double visc_valy = 0.0;
      double visc_valz = 0.0;
      double vms1 = 0.0;
      double vms2 = 0.0;
      double vms3 = 0.0;
      double supg_term = 0.0;
      double lsic_term = 0.0;
      gradu.Mult(rm,gradu_rm); // gradu * rm
      for (int a = 0; a < nd_u; ++a)
      {
         dshapedx_u.GetRow(a,Nadx);
         if(viscous_term == ifs_solver::ViscousTermType::LaplaceForm)
         {
            gradu.Mult(Nadx,Nadx_sigma); // Nadx_sigma = gradv_ij * N^a,j
            visc_valx = -p * Nadx(0) + mu_f * Nadx_sigma(0);
            visc_valy = -p * Nadx(1) + mu_f * Nadx_sigma(1);
            if(dim == 3) visc_valz = -p * Nadx(2) + mu_f * Nadx_sigma(2);
         }
         else if (viscous_term == ifs_solver::ViscousTermType::StressForm)
         {
            fluid_stress.Mult(Nadx,Nadx_sigma);
            visc_valx = Nadx_sigma(0);
            visc_valy = Nadx_sigma(1);
            if(dim == 3) visc_valz = Nadx_sigma(2);
         }
         else
         {
            MFEM_ABORT("Viscous term type not supported");
         }

         /////// STABLE METHOD ///////
         (*elvec[0])(a)                        += w * ( rm_inertial(0)*shape_u(a) + visc_valx);
         (*elvec[0])(a + nd_u)                 += w * ( rm_inertial(1)*shape_u(a) + visc_valy);
         if (dim == 3) (*elvec[0])(a + 2*nd_u) += w * ( rm_inertial(2)*shape_u(a) + visc_valz);

         /////// STABILIZATIONS ///////
         if (supg_flag)
         {
            // supg term
            supg_term = (Nadx*u) * w * tau_m;

            (*elvec[0])(a)                        += supg_term * rm(0);
            (*elvec[0])(a + nd_u)                 += supg_term * rm(1);
            if (dim == 3) (*elvec[0])(a + 2*nd_u) += supg_term * rm(2);
         }
         if (vms_flag)
         {
            // vms term 1
            vms1 = w * -shape_u(a) * tau_m;

            // vms term 2
            vms2 = w * tau_m * tau_m / rho * -(Nadx * rm);

            // additional momentum restoring term
            vms3 = w * tau_m * tau_m * tau_b * (Nadx * rm);

            (*elvec[0])(a)                        += vms2 * rm(0) + (vms1 + vms3) * gradu_rm(0);
            (*elvec[0])(a + nd_u)                 += vms2 * rm(1) + (vms1 + vms3) * gradu_rm(1);
            if (dim == 3) (*elvec[0])(a + 2*nd_u) += vms2 * rm(2) + (vms1 + vms3) * gradu_rm(2);
         }
         if (lsic_flag)
         {
            // lsic term
            lsic_term = w * tau_c * rho * divu;
            
            (*elvec[0])(a)                        += lsic_term * Nadx(0);
            (*elvec[0])(a + nd_u)                 += lsic_term * Nadx(1);
            if (dim == 3) (*elvec[0])(a + 2*nd_u) += lsic_term * Nadx(2);
         }
      }

      ///////// continuity residual /////////
      elvec[1]->Add(w*divu,shape_p);     
      if (pspg_flag)
      {
         dshapedx_p.AddMult(rm,*elvec[1],w * tau_m / rho);
      }
   }
}

void ifs_solver::IncompressibleFluidIntegrator::AssembleElementGrad(const Array<const FiniteElement*> &el,
                                     ElementTransformation &trans,
                                     const Array<const Vector *> &elfun,
                                     const Array2D<DenseMatrix *> &elmats)
{
   if (el.Size() != 2)
   {
      mfem_error("Finite element space of incorrect block number");
   }

   ///////// setup data /////////
   int nd_u = el[0]->GetDof();
   int nd_p = el[1]->GetDof();

   int dim = el[0]->GetDim();
   int spaceDim = trans.GetSpaceDim();

   if (dim != spaceDim)
   {
      mfem_error("Not defined on manifold meshes");
   }

   elmats(0,0)->SetSize(nd_u*dim, nd_u*dim);
   elmats(0,1)->SetSize(nd_u*dim, nd_p);
   elmats(1,0)->SetSize(nd_p, nd_u*dim);
   elmats(1,1)->SetSize(nd_p, nd_p);

   *elmats(0,0) = 0.0;
   *elmats(0,1) = 0.0;
   *elmats(1,0) = 0.0;
   *elmats(1,1) = 0.0;

   dshape_u.SetSize(nd_u, dim);
   dshapedx_u.SetSize(nd_u, dim);
   shape_u.SetSize(nd_u);
   gradu.SetSize(dim);
   u.SetSize(dim);
   u_dot_gradu.SetSize(dim);
   EF.UseExternalData(elfun[0]->GetData(),nd_u,dim); // velocity dofs
   shape_u_lap.SetSize(nd_u);
   shape_u_lap = 0.0;
   Nadx_sigma.SetSize(dim);
   u_lap.SetSize(dim);
   u_lap = 0.0;
   if (mesh_velo) u_mesh.SetSize(dim);

   dshape_p.SetSize(nd_p, dim);
   shape_p.SetSize(nd_p);
   dshapedx_p.SetSize(nd_p, dim);
   gradp.SetSize(dim);
   
   Nadx.SetSize(dim);
   Nbdx.SetSize(dim);

   rm.SetSize(dim);
   rm_inertial.SetSize(dim);

   fluid_stress.SetSize(dim);
   f_vec.SetSize(dim);
   
   double w = 0.0;
   double rho = 0.0;

   // time factors
   double am = alpha_m;
   double af = alpha_f * gamma * (*dt);

   ///////// integration rule /////////
   const IntegrationRule *ir = nullptr;
   if (ir == nullptr)
   {
      int order = 2*el[0]->GetOrder(); // exact integration for linear terms
      ir = &IntRules.Get(el[0]->GetGeomType(), order);
   }

   ///////// loop over integration points /////////
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      ///////// compute data /////////
      const IntegrationPoint &ip = ir->IntPoint(i);
      trans.SetIntPoint(&ip);

      // pressure shape functions
      el[1]->CalcShape(ip,shape_p);
      el[1]->CalcDShape(ip, dshape_p);
      Mult(dshape_p, trans.InverseJacobian(), dshapedx_p);
      dshapedx_p.MultTranspose(*elfun[1],gradp); // grad p (check this!)

      // velocity shape functions
      el[0]->CalcShape(ip,shape_u);
      el[0]->CalcDShape(ip,dshape_u);
      Mult(dshape_u, trans.InverseJacobian(), dshapedx_u);
      MultAtB(EF, dshapedx_u, gradu); // grad u
      EF.MultTranspose(shape_u,u); // u
      if (mesh_velo)
      {
         mesh_velo->Eval(u_mesh,trans,ip);
         u -= u_mesh; // u = (u - u_mesh)
      }
      gradu.Mult(u,u_dot_gradu); // u dot gradu

      if (nxx_flag)
      {
         // Not needed for p=1, minimal effect if included
         // this will be zero on linear triangles and tets
         // H1FECollection: Will work for all element types except wedges and pyramids
         el[0]->CalcPhysLaplacian(trans,shape_u_lap);
      }      

      double mu_f = dyn_vis->Eval(trans,ip);
      if (vms_flag)
      {
         EF.MultTranspose(shape_u_lap,u_lap); // u_i,jj

         ///// Compute rm at the integration point /////
         // rm = rho * (udot + u*gradu - f) + gradp - mu*d2u/dxdx
         Compute_rm(rm,rm_inertial,mu_f,trans,ip,el);
      }      

      // stabilization parameters
      if (supg_flag) tau_m = tau_m_evaluator->EvalTauM(trans,ip,u);
      if (vms_flag) tau_b = tau_m_evaluator->EvalTauB(trans,ip,tau_m,rm);
      if (lsic_flag) tau_c = mass_stab_coeff->Eval_Coeff(tau_m); // o.w tau_c = 0.0

      ///// Add contributions to element vector /////
      w = ip.weight * trans.Weight(); // w = w_q * detJ
      double rho = rhof->Eval(trans,ip);

      ///// (0,0) block /////
      double Nadx_u = 0.0;
      double Nbdx_u = 0.0;
      double Nbdx_rm = 0.0;
      double Nadx_rm = 0.0;

      double ri = 0.0;
      double rv_1 = 0.0;
      double rv_2 = 0.0;

      double supg_term = 0.0;
      double vms1 = 0.0;
      double vms2 = 0.0;
      double vms3 = 0.0;
      double lsic_term = 0.0;

      for (int a=0; a < nd_u; ++a)
      {
         dshapedx_u.GetRow(a,Nadx);
         Nadx_u = Nadx * u;
         Nadx_rm = Nadx * rm;

         for (int b=0; b < nd_u; ++b)
         {
            dshapedx_u.GetRow(b,Nbdx);
            Nbdx_u = Nbdx * u;
            Nbdx_rm = Nbdx * rm;

            // rho * (vdot + v*gradv)
            ri = w*rho*shape_u(a)*(am * shape_u(b) + af*Nbdx_u);

            // Nax*Nbx
            if(viscous_term == ifs_solver::ViscousTermType::LaplaceForm)
            {
               rv_1 = w*af*mu_f*(Nadx*Nbdx);
               rv_2 = 0.0;
            }
            else if (viscous_term == ifs_solver::ViscousTermType::StressForm)
            {
               rv_1 = w*af*mu_f*(Nadx*Nbdx);
               rv_2 = w*af*mu_f;
            }
            else
            {
               MFEM_ABORT("Viscous term type not supported");
            }

            /////// STABLE METHOD ///////
            (*elmats(0,0))(a,b)                               += ri + rv_1 + rv_2 * Nadx(0) * Nbdx(0);
            (*elmats(0,0))(a,b + nd_u)                        += rv_2 * Nadx(1) * Nbdx(0);
            if (dim == 3) (*elmats(0,0))(a,b + 2*nd_u)        += rv_2 * Nadx(2) * Nbdx(0);
            
            (*elmats(0,0))(a + nd_u,b)                        += rv_2 * Nadx(0) * Nbdx(1);
            (*elmats(0,0))(a + nd_u, b + nd_u)                += ri + rv_1 + rv_2 * Nadx(1) * Nbdx(1);
            if (dim == 3) (*elmats(0,0))(a + nd_u,b + 2*nd_u) += rv_2 * Nadx(2) * Nbdx(1);
            
            if (dim == 3)
            {
               (*elmats(0,0))(a + 2*nd_u,b)                   += rv_2 * Nadx(0) * Nbdx(2);
               (*elmats(0,0))(a + 2*nd_u,b + nd_u)            += rv_2 * Nadx(1) * Nbdx(2);
               (*elmats(0,0))(a + 2*nd_u, b + 2*nd_u)         += ri + rv_1 + rv_2 * Nadx(2) * Nbdx(2);
            }

            /////// STABILIZATIONS ///////
            if (supg_flag)
            {
               // supg term
               supg_term = w * tau_m * Nadx_u * (rho * ( am * shape_u(b) + Nbdx_u * af ) - af * mu_f * shape_u_lap(b));

               (*elmats(0,0))(a,b)                                  += supg_term;
               (*elmats(0,0))(a + nd_u, b + nd_u)                   += supg_term;
               if (dim == 3) (*elmats(0,0))(a + 2*nd_u, b + 2*nd_u) += supg_term;
            }
            if (vms_flag)
            {
               // vms term 1
               vms1 = w * af * tau_m * Nbdx_rm * -shape_u(a);

               // vms term 2
               vms2 = -w * tau_m * tau_m * Nadx_rm * (am * shape_u(b) + af * Nbdx_u - af * mu_f * shape_u_lap(b) / rho); 

               // additional momentum restoring term
               vms3 = w * af * tau_m * Nadx_rm * tau_b * Nbdx_rm * tau_m;

               (*elmats(0,0))(a,b)                                  += vms1 + vms2 + vms3;
               (*elmats(0,0))(a + nd_u, b + nd_u)                   += vms1 + vms2 + vms3;
               if (dim == 3) (*elmats(0,0))(a + 2*nd_u, b + 2*nd_u) += vms1 + vms2 + vms3;
            }
            if (lsic_flag)
            {
               // lsic term
               lsic_term = af * w * tau_c * rho;

               (*elmats(0,0))(a,b)                               += lsic_term * Nadx(0) * Nbdx(0);
               (*elmats(0,0))(a,b + nd_u)                        += lsic_term * Nadx(0) * Nbdx(1);
               if (dim == 3) (*elmats(0,0))(a,b + 2*nd_u)        += lsic_term * Nadx(0) * Nbdx(2);
               
               (*elmats(0,0))(a + nd_u,b)                        += lsic_term * Nadx(1) * Nbdx(0);
               (*elmats(0,0))(a + nd_u, b + nd_u)                += lsic_term * Nadx(1) * Nbdx(1);
               if (dim == 3) (*elmats(0,0))(a + nd_u,b + 2*nd_u) += lsic_term * Nadx(1) * Nbdx(2);
               
               if (dim == 3)
               {
                  (*elmats(0,0))(a + 2*nd_u,b)                   += lsic_term * Nadx(2) * Nbdx(0);
                  (*elmats(0,0))(a + 2*nd_u,b + nd_u)            += lsic_term * Nadx(2) * Nbdx(1);
                  (*elmats(0,0))(a + 2*nd_u, b + 2*nd_u)         += lsic_term * Nadx(2) * Nbdx(2);
               }

            }
         }
      }  

      ///// (0,1) and (1,0) block /////
      double Np_w_af = 0.0;
      double pspg_term = 0.0;
      for (int a = 0; a < nd_u; ++a)
      {
         dshapedx_u.GetRow(a,Nadx);
         Nadx_u = Nadx * u;
         vms2 = (Nadx * rm) * tau_m / rho;

         for (int b = 0; b < nd_p; ++b)
         {
            ///// (0,1) block /////
            Np_w_af = -w * af * shape_p(b);

            /////// STABLE METHOD ///////
            (*elmats(0,1))(a,b)                        += Np_w_af * Nadx(0);
            (*elmats(0,1))(a + nd_u,b)                 += Np_w_af * Nadx(1);
            if (dim == 3) (*elmats(0,1))(a + 2*nd_u,b) += Np_w_af * Nadx(2);

            /////// STABILIZATIONS ///////
            if (supg_flag || vms_flag)
            {
               supg_term = supg_flag ? w * af * tau_m * Nadx_u : 0.0;
               vms2 = vms_flag ? -w * tau_m * af : 0.0; 

               (*elmats(0,1))(a,b)                        += (supg_term + vms2) * dshapedx_p(b,0);
               (*elmats(0,1))(a + nd_u,b)                 += (supg_term + vms2) * dshapedx_p(b,1);
               if (dim == 3) (*elmats(0,1))(a + 2*nd_u,b) += (supg_term + vms2) * dshapedx_p(b,2);
            }

            ///// (1,0) block /////
            Np_w_af *= -1.0;

            /////// STABLE METHOD ///////
            (*elmats(1,0))(b,a)                        += Np_w_af * Nadx(0);
            (*elmats(1,0))(b,a + nd_u)                 += Np_w_af * Nadx(1);
            if (dim == 3) (*elmats(1,0))(b,a + 2*nd_u) += Np_w_af * Nadx(2);

            /////// STABILIZATIONS ///////
            if (pspg_flag)
            {
               pspg_term = w * tau_m * (am * shape_u(a) + af * Nadx_u - af * mu_f * shape_u_lap(a) / rho );

               (*elmats(1,0))(b,a)                        += pspg_term * dshapedx_p(b,0);
               (*elmats(1,0))(b,a + nd_u)                 += pspg_term * dshapedx_p(b,1);
               if (dim == 3) (*elmats(1,0))(b,a + 2*nd_u) += pspg_term * dshapedx_p(b,2);
            }
         }
      }

      ///// (1,1) block /////
      if (pspg_flag)
      {
         AddMult_a_AAt(w * tau_m * af / rho,dshapedx_p,*elmats(1,1));
      }
   }
}

void ifs_solver::IncompressibleFluidIntegrator::Compute_rm(Vector &rm, Vector &rm_inertial, double &mu_f, ElementTransformation &trans, const IntegrationPoint &ip, const Array<const FiniteElement *> &el)
{
   ///// compute Rm at integration point /////

   // time term
   vdot_gf->GetVectorValue(trans,ip,rm);

   // nonlinear term
   rm.Add(1.0,u_dot_gradu);

   // forcing (assumes forcing is applied over entire domain, which is assumed in general)
   for (auto &bfs : *body_force_terms)
   { 
      bfs.coeff->Eval(f_vec,trans,ip);
      rm.Add(-1.0,f_vec);
   }

   // mult by density
   if (rhof)
   {
      rm *= rhof->Eval(trans,ip);
   }

   rm_inertial.Set(1.0,rm);

   // grad p
   rm.Add(1.0,gradp);   

   // viscous term
   rm.Add(-1.0 * mu_f,u_lap);
}