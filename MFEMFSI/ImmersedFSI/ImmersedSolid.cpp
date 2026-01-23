#include "ImmersedSolid.hpp"
#include "NodalAverageOutputs.h"
#include <general/forall.hpp>
#include "../Utilities/linear_element_ext.hpp"
#ifdef MFEM_USE_MPI
#include "RefBoundaryTractionIntegrator.hpp"
#include "ImmersedFollowerPressureLoad.hpp"
#endif

using namespace mfem;
using namespace std;
using namespace mat;
using namespace immersed_solid;

#ifdef MFEM_USE_MPI
ImmersedSolid::ImmersedSolid(ParMesh *solid_mesh, int order_, Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
   ParFiniteElementSpace *vfes_, ParMesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_)
#else
ImmersedSolid::ImmersedSolid(Mesh *solid_mesh, int order_, Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
   FiniteElementSpace *vfes_, Mesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_)
#endif
: pmesh(solid_mesh), orderd(order_), init_rhos(init_rho_s_), fluid_density(fluid_density_), mu_f(mu_f_), vfes(vfes_), fluid_pmesh(fluid_mesh), 
mat_type(mat_), orderingtype(orderingtype_), formulation_type(form_type_), subtract_fluid_visc(subtract_visc_), accel_terms(accel_)
{
   int dim = pmesh->Dimension();
   // finite element collection (Lagrange polynomials thru Gauss-Lobatto points)
   dfec = new H1_FECollection(orderd,dim);

   // finite element space
   #ifdef MFEM_USE_MPI
   dfes = new ParFiniteElementSpace(pmesh,dfec,dim,orderingtype);
   #else
   dfes = new FiniteElementSpace(pmesh, dfec, dim, orderingtype);
   #endif

   // set mesh nodal gridfunction
   pmesh->SetNodalFESpace(dfes);   
   
   // setup boundary attribute arrays
   // displacement bcs
   disp_attr.SetSize(pmesh->bdr_attributes.Max()*dim);
   disp_attr = 0;

   #ifdef MFEM_USE_MPI
   // traction bcs
   traction_attr.SetSize(pmesh->bdr_attributes.Max());
   traction_attr = 0;
   #endif

   int dfes_truevsize = dfes->GetTrueVSize();

   // setup and initialize variables
   dn_gf.SetSpace(dfes);
   dn_gf = 0.0;
   ddotn_gf.SetSpace(dfes);
   ddotn_gf = 0.0;
   d2dotn_gf.SetSpace(dfes);
   d2dotn_gf = 0.0;
   ddotn.SetSize(dfes_truevsize);
   ddotn = 0.0;
   ddotnp1.SetSize(dfes_truevsize);
   ddotnp1 = 0.0;
   dg.SetSize(dfes_truevsize);
   dg = 0.0;
   ddotg.SetSize(dfes_truevsize);
   ddotg = 0.0;
   dn.SetSize(dfes_truevsize);
   dn = 0.0;
   dnp1.SetSize(dfes_truevsize);
   dnp1 = 0.0;   

   temp_del_gf.SetSpace(vfes);
   temp_del_gf = 0.0;
   temp_delvdot_ib_gf.SetSpace(dfes);
   temp_delvdot_ib_gf = 0.0;
   delvdot_ib.SetSize(dfes_truevsize);
   delvdot_ib = 0.0;
   dtemp_gf.SetSpace(dfes);
   dtemp_gf = 0.0;

   // Store initial mesh positions.
   orig_nodes.SetSpace(dfes);
   #ifdef MFEM_USE_MPI
   orig_nodes = *dynamic_cast<ParGridFunction *>(pmesh->GetNodes());
   #else
   orig_nodes = *pmesh->GetNodes();
   #endif
   new_nodes.SetSpace(dfes);
   new_nodes = 0.0;

   d2dotn.SetSize(dfes_truevsize);
   d2dotn = 0.0;
   d2dotnp1.SetSize(dfes_truevsize);
   d2dotnp1 = 0.0;
   d2dotg.SetSize(dfes_truevsize);
   d2dotg = 0.0;
   density_diff = new SumCoefficient(fluid_density,*init_rhos,-1.0,1.0);
   res.SetSize(dfes_truevsize);
   res = 0.0;
   ft.SetSize(dfes_truevsize);
   ft = 0.0;

   #ifdef MFEM_USE_MPI
   my_rank = vfes->GetMyRank();
   #else
   my_rank = 0;
   #endif

   // setup point search
   #ifdef MFEM_USE_GSLIB   

   #ifdef MFEM_USE_MPI
   finder_fluid = new ImmersedFSIFindPts(vfes->GetComm());
   #else
   finder_fluid = new ImmersedFSIFindPts();
   #endif

   finder_fluid->Setup(*fluid_pmesh);

   #ifdef MFEM_USE_MPI
   finder_fluid->SetupLocalDofToMasterRankMap(vfes);
   finder_fluid->SetupLocalTDofToLDofSolidMap(dfes);
   #endif
   #endif

   true_interp_points_in_solid.SetSize(dfes->GetTrueVSize());
   true_interp_points_in_solid = 0.0;

   interp_points_in_solid_tn.SetSize(dfes->GetVSize());
   interp_points_in_solid_tn = 0.0;

   interp_points_in_solid = *pmesh->GetNodes();
   GridFunction *nodes = pmesh->GetNodes();
   #ifdef MFEM_USE_MPI
   dynamic_cast<ParGridFunction *>(nodes)->ParallelProject(true_interp_points_in_solid);
   #else
   nodes->GetTrueDofs(true_interp_points_in_solid);
   #endif
}

void ImmersedSolid::StabsScaling(mfem::GridFunction *s_gf_fluid_)
{
   stabs_scaling = true;
   s_gf_fluid = s_gf_fluid_;

   // setup VMS stab scaling
   if (stabs_scaling)
   {
      // setup interpolation points (only boundary dofs of the solid)
      Array<int> bdr_attr(pmesh->bdr_attributes.Max());
      bdr_attr = 1;
      dfes->GetEssentialTrueDofs(bdr_attr,solid_bdr_tdofs);
      solid_bdr_nodes.SetSize(solid_bdr_tdofs.Size());
      solid_bdr_nodes = 0.0;
      solid_bdr_nodes_tn.SetSize(solid_bdr_tdofs.Size());
      solid_bdr_nodes_tn = 0.0;

      const int Ndofs = solid_bdr_tdofs.Size();
      const auto id = solid_bdr_tdofs.Read();
      auto FV_nodes = true_interp_points_in_solid.Read();
      auto FV = solid_bdr_nodes.ReadWrite();
      mfem::forall(Ndofs, [=] MFEM_HOST_DEVICE (int i) { FV[i] = FV_nodes[id[i]]; });
   }
}

void ImmersedSolid::Setup(double dt, double alpha_m_, double alpha_f_, double gamma_, double beta_)
{
   if (my_rank == 0)
   {
      mfem::out << "Setup Immersed Solid" << std::endl;     
   }
   PrintInfo();

   curr_d2dot.SetSpace(dfes);
   curr_d2dot = 0.0;
   curr_ddot.SetSpace(dfes);
   curr_ddot = 0.0;

   #ifdef MFEM_USE_MPI
   iSolidNLF = new ParNonlinearForm(dfes);
   #else
   iSolidNLF = new NonlinearForm(dfes);
   #endif
   this->dt = dt;
   iSolidNLF->AddDomainIntegrator(new ImmersedSolidIntegrator(mat_type,formulation_type,&orig_nodes,&curr_d2dot,&curr_ddot,&this->dt,beta_,gamma_,alpha_m_,alpha_f_,accel_terms,init_rhos,fluid_density,mu_f,subtract_fluid_visc));      
   
   // BCs
   #ifdef MFEM_USE_MPI
   // follower pressure bcs
   minusOne.constant = -1.0; // to account for sign difference in follower load vs stress term in eqns
   for (auto &fp_terms : fp_bcs)
   {
      iSolidNLF->AddBoundaryIntegrator(new ImmersedFollowerPressureLoad(formulation_type,fp_terms.coeff,&minusOne,&this->dt,alpha_f_,beta_),fp_terms.attr);      
   }

   // traction bcs (no minus sign)
   Ft_form = new ParLinearForm(dfes);
   for (auto &trac_bc : traction_bcs)
   {
      Ft_form->AddBoundaryIntegrator(new RefBoundaryTractionIntegrator(*trac_bc.coeff,formulation_type,&orig_nodes),trac_bc.attr);
   }
   #endif

   beta = beta_;
   gamma = gamma_;
   alpha_m = alpha_m_;
   alpha_f = alpha_f_;

   // set initial conditions
   #ifdef MFEM_USE_MPI
   dn_gf.ParallelProject(dn);
   #else
   dn_gf.GetTrueDofs(dn);
   #endif
}

void ImmersedSolid::SetTimeAndDt(double time_, double dt_)
{
   time = time_;
   dt = dt_;
}

void ImmersedSolid::SetWeakBCs(const double &time)
{
   #ifdef MFEM_USE_MPI
   // set time on traction bcs
   for (auto &trac_bc : traction_bcs)
   {
      trac_bc.coeff->SetTime(time);
   }
   Ft_form->Assemble();
   Ft_form->ParallelAssemble(ft); // evaluate at new time t
   #endif

   // update body force fb, this is handled via flow solver
   // set time on body force bcs
   for (auto &bf_term : *accel_terms)
   {
      bf_term.coeff->SetTime(time);
   }     

   // update follower pressure loads
   for (auto &fp_terms : fp_bcs)
   {
      fp_terms.coeff->SetTime(time);
   }
}

void ImmersedSolid::ComputeSolidResidualAndSpread2Fluid(Vector &lambda_fsi_vec)
{
   // Compute solid residual at midpoint configuration
   // Rs = M*anp1 + Fint - Fext

   sw_res->Start();
   // set solid velocity and acceleration grid function
   curr_d2dot.SetFromTrueDofs(d2dotg);
   curr_ddot.SetFromTrueDofs(ddotg);

   // compute immersed solid residual
   iSolidNLF->Mult(dg,res); 
   
   #ifdef MFEM_USE_MPI
   // traction bcs
   res.Add(-1.0, ft);
   #endif
   sw_res->Stop();

   sw_fsires->Start();
   #ifdef MFEM_USE_GSLIB
   #ifdef MFEM_USE_MPI
   HypreParMatrix *C = dynamic_cast<HypreParMatrix *>(&finder_fluid->GetCouplingMatrix());
   #else
   SparseMatrix *C = dynamic_cast<SparseMatrix *>(&finder_fluid->GetCouplingMatrix());
   #endif
   #else
   #ifdef MFEM_USE_MPI
   HypreParMatrix *C = nullptr;
   #else
   SparseMatrix *C = nullptr;
   #endif
   #endif

   C->Mult(res,lambda_fsi_vec);

   // #ifdef MFEM_USE_GSLIB
   // finder_fluid->SpreadForces(true_interp_points_in_solid,res,temp_del_gf,orderingtype);
   // #ifdef MFEM_USE_MPI
   // temp_del_gf.ParallelProject(lambda_fsi_vec);
   // #else
   // temp_del_gf.GetTrueDofs(lambda_fsi_vec);
   // #endif
   // #endif

   sw_fsires->Stop();
}

Operator &ImmersedSolid::FormLHS()
{
   sw_tan->Start();
   // solid tangent matrix
   #ifdef MFEM_USE_MPI
   HypreParMatrix *Kmat = dynamic_cast<HypreParMatrix *>(&iSolidNLF->GetGradient(dg)); 
   #else
   SparseMatrix *Kmat = dynamic_cast<SparseMatrix *>(&iSolidNLF->GetGradient(dg)); 
   #endif

   sw_tan->Stop();

   sw_fsitan->Start();
   #ifdef MFEM_USE_MPI
   // Get C matrix (fluid x solid dofs)
   #ifdef MFEM_USE_GSLIB
   HypreParMatrix *C = dynamic_cast<HypreParMatrix *>(&finder_fluid->GetCouplingMatrix());
   HypreParMatrix *Ct = dynamic_cast<HypreParMatrix *>(&finder_fluid->GetCouplingMatrixTranspose());
   #else
   HypreParMatrix *C = nullptr;
   HypreParMatrix *Ct = nullptr;
   #endif

   // C K Ct
   if (KCt) delete KCt;
   KCt = ParMult(Kmat,Ct,true);
   if (CKCt) delete CKCt;
   CKCt = ParMult(C,KCt,true);

   sw_fsitan->Stop();
   return *CKCt;
   
   #else
   // Get C matrix (fluid x solid dofs)
   #ifdef MFEM_USE_GSLIB
   SparseMatrix *C = dynamic_cast<SparseMatrix *>(&finder_fluid->GetCouplingMatrix());
   SparseMatrix *Ct = dynamic_cast<SparseMatrix *>(&finder_fluid->GetCouplingMatrixTranspose());
   #else
   SparseMatrix *C = nullptr;
   SparseMatrix *Ct = nullptr;
   #endif

   // C K Ct
   if (KCt) delete KCt;
   KCt = mfem::Mult(*Kmat,*Ct);
   if (CKCt) delete CKCt;
   CKCt = mfem::Mult(*C,*KCt);

   sw_fsitan->Stop();
   return *CKCt;
   #endif
}

void ImmersedSolid::Predictors(Vector &vf, Vector &vdotf)
{
   // set ddotnp1 and d2dotnp1 from fluid predictors for consistency!
   // ddotnp1 = v_np1
   temp_del_gf.SetFromTrueDofs(vf);
   #ifdef MFEM_USE_GSLIB
   finder_fluid->Interpolate(interp_points_in_solid,temp_del_gf,temp_delvdot_ib_gf,orderingtype);
   #endif
   #ifdef MFEM_USE_MPI
   temp_delvdot_ib_gf.ParallelProject(ddotnp1);
   #else
   temp_delvdot_ib_gf.GetTrueDofs(ddotnp1);
   #endif

   // d2dotnp1 = vdot_np1
   temp_del_gf.SetFromTrueDofs(vdotf);
   #ifdef MFEM_USE_GSLIB
   finder_fluid->Interpolate(interp_points_in_solid,temp_del_gf,temp_delvdot_ib_gf,orderingtype);
   #endif
   #ifdef MFEM_USE_MPI
   temp_delvdot_ib_gf.ParallelProject(d2dotnp1);
   #else
   temp_delvdot_ib_gf.GetTrueDofs(d2dotnp1);   
   #endif

   auto anp1_R = d2dotnp1.Read();
   auto an_R = d2dotn.Read();
   auto vn_R = ddotn.Read(); 
   auto dnp1_RW = dnp1.ReadWrite();
   auto dn_R = dn.Read();
   const int N = dnp1.Size();
   const double dt_ = dt;
   const double fac_n = 0.5*dt*dt*(1.0 - 2.0*beta);
   const double fac_np1 = 0.5*dt*dt*(2.0*beta);
   mfem::forall(N, [=] MFEM_HOST_DEVICE (int i) 
   { 
      dnp1_RW[i] = dn_R[i] + dt_*vn_R[i] + fac_n*an_R[i] + fac_np1*anp1_R[i]; 
   });

   // store tn position of solid
   interp_points_in_solid_tn = interp_points_in_solid;
   if (stabs_scaling) solid_bdr_nodes_tn = solid_bdr_nodes;
}

void ImmersedSolid::SetEssentialBCs(const double &time)
{
   // nothing to do
}

void ImmersedSolid::Initiator()
{
   add((1.0 - alpha_f),dn,alpha_f,dnp1,dg);
   add((1.0 - alpha_f),ddotn,alpha_f,ddotnp1,ddotg);
   add((1.0 - alpha_m),d2dotn,alpha_m,d2dotnp1,d2dotg);

   // move mesh and reset point search
   dtemp_gf.SetFromTrueDofs(dg);
   MoveMesh(dtemp_gf); // contains current displacements
   // UpdatePointSearch(dtemp_gf,dg); // not needed since motion map for lambda terms fixed at tn

   if (stabs_scaling)
   {
      #ifdef MFEM_USE_GSLIB
      MFEM_ASSERT(dfes->GetOrdering() == s_gf_fluid->FESpace()->GetOrdering(),"Ordering of solid and fluid gridfunctions must match for special VMS stabs scaling");
      finder_fluid->TransferIndicator(solid_bdr_nodes_tn,*s_gf_fluid,dfes->GetOrdering());
      // note internally s_gf_fluid is zeroed prior to interpolation of s_gf to its fluid counterpart
      #endif
   }
}

void ImmersedSolid::Correctors(Vector &del)
{
   // project delvdot to solid mesh
   temp_del_gf.SetFromTrueDofs(del);
   #ifdef MFEM_USE_GSLIB
   finder_fluid->Interpolate(interp_points_in_solid_tn,temp_del_gf,temp_delvdot_ib_gf,orderingtype);
   #endif
   #ifdef MFEM_USE_MPI
   temp_delvdot_ib_gf.ParallelProject(delvdot_ib);
   #else
   temp_delvdot_ib_gf.GetTrueDofs(delvdot_ib);
   #endif
   
   d2dotnp1.Add(-1.0,delvdot_ib);
   ddotnp1.Add(-gamma*dt,delvdot_ib);
   dnp1.Add(-beta*dt*dt,delvdot_ib);
}

void ImmersedSolid::Update()
{
   dn = dnp1;
   dn_gf.SetFromTrueDofs(dn);

   ddotn = ddotnp1;
   ddotn_gf.SetFromTrueDofs(ddotn);

   d2dotn = d2dotnp1;
   d2dotn_gf.SetFromTrueDofs(d2dotn);

   // update solid mesh (if UL formulation) //
   MoveMesh(dn_gf);
   UpdatePointSearch(dn_gf,dn);
   if (stabs_scaling)
   {
      #ifdef MFEM_USE_GSLIB
      finder_fluid->TransferIndicator(solid_bdr_nodes,*s_gf_fluid,dfes->GetOrdering());
      #endif
   }
}

void ImmersedSolid::MoveMesh(GridFunction &total_disp_gf)
{
   // update solid mesh (if UL formulation) //
   if (formulation_type == form_type::UpdatedLagrangian)
   {
      // // move mesh
      // new_nodes = orig_nodes; // x0
      // new_nodes += dn_gf; // x_new = disp + x0, where disp is disp relative to ref config
      add(1.0, orig_nodes, 1.0, total_disp_gf, new_nodes);

      // set mesh to new nodes
      pmesh->NewNodes(new_nodes);   
   }
}

void ImmersedSolid::UpdatePointSearch(const GridFunction &d, const Vector &true_d)
{
   /// Update solid points ///
   // update interpolation nodes from fluid->solid
   interp_points_in_solid = *pmesh->GetNodes();
   if (formulation_type == form_type::TotalLagrangian)
   {
      interp_points_in_solid.Add(1.0,d);
   }

   #ifdef MFEM_USE_MPI
   ParGridFunction *nodes = dynamic_cast<ParGridFunction *>(pmesh->GetNodes());
   nodes->ParallelProject(true_interp_points_in_solid);
   #else
   GridFunction *nodes = dynamic_cast<GridFunction *>(pmesh->GetNodes());
   nodes->GetTrueDofs(true_interp_points_in_solid);
   #endif
   if (formulation_type == form_type::TotalLagrangian)
   {
      true_interp_points_in_solid.Add(1.0,true_d);
   }

   if (stabs_scaling)
   {
      const int Ndofs= solid_bdr_tdofs.Size();
      const auto id = solid_bdr_tdofs.Read();
      auto FV_nodes = true_interp_points_in_solid.Read();
      auto FV = solid_bdr_nodes.ReadWrite();
      mfem::forall(Ndofs, [=] MFEM_HOST_DEVICE (int i) { FV[i] = FV_nodes[id[i]]; });
   }

   /// Reset gslib class data ///
   #ifdef MFEM_USE_GSLIB
   finder_fluid->Setup(*fluid_pmesh);
   #endif
}

void ImmersedSolid::TrackNodalDisplacement(Vector &coord, Vector &disp, double tol)
{
   if (find_node)
   {
      // find the node number with position provided by coords
      #ifdef MFEM_USE_GSLIB
      #ifdef MFEM_USE_MPI
      FindPointsGSLIB finder(dfes->GetComm());
      #else
      FindPtsGSLIB finder;
      #endif
      finder.FindPoints(*pmesh,coord);
      Array<int> code = finder.GetCode();
      Array<int> elem_ids_gslib = finder.GetElem();
      Array<int> proc_ids_gslib = finder.GetProc();
      MFEM_ASSERT(code[0] != 2, "Error: node to track was not found!");
      int elem_num = elem_ids_gslib[0];
      int proc_id = proc_ids_gslib[0];
      #else
      DenseMatrix point_mat(pmesh->Dimension(),1);
      point_mat.SetCol(0,coord);
      point_mat.Print();
      Array<int> elem_ids;
      Array<IntegrationPoint> ips;
      int npts_found = pmesh->FindPoints(point_mat,elem_ids,ips);
      MFEM_ASSERT(npts_found == 1, "Error: node to track was not found!");
      int elem_num = elem_ids[0];
      #endif

      local_node_num = -1;
      #ifdef MFEM_USE_GSLIB
      if (proc_id == my_rank)
      {
      #endif
         int dim = pmesh->Dimension();
         Vector curr_ip(dim);
         if (elem_num > -1)
         {
            const FiniteElement *pfe = dfes->GetFE(elem_num);
            ElementTransformation *eltrans = dfes->GetElementTransformation(elem_num);
            Array<int> dofs;
            auto doftrans = dfes->GetElementDofs(elem_num,dofs);
            // integration rule = nodes of element
            // i.e. int point 0 = node 0
            const IntegrationRule *ir = &pfe->GetNodes();

            double diff_x = 1.0;
            double diff_y = 1.0;
            double diff_z = 1.0;

            for (int j = 0; j < ir->GetNPoints(); j++)
            {
               const IntegrationPoint &ip = ir->IntPoint(j);

               eltrans->Transform(ip,curr_ip);

               // compare nodes with coords
               diff_x = abs(curr_ip(0) - coord(0));
               diff_y = abs(curr_ip(1) - coord(1));
               if (dim == 3) diff_z = abs(curr_ip(2) - coord(2));

               if (diff_x <= tol && diff_y <= tol)
               {
                  if (dim == 3)
                  {
                     if (diff_z <= tol)
                     {
                        local_node_num = dofs[j];
                        break;
                     }
                  }
                  else
                  {
                     local_node_num = dofs[j];
                     break;
                  }
               }
            }

            MFEM_VERIFY(local_node_num > -1, "Error: node not found!");
         }
      #ifdef MFEM_USE_GSLIB
      }
      #endif
      find_node = false;
   }

   // return displacement
   int id0, id1, id2;
   double x,y,z,x_glob,y_glob,z_glob;
   x = y = z = 0.0;
   x_glob = y_glob = z_glob = 0.0;
   int dim = pmesh->Dimension();
   if (local_node_num > -1)
   {
      id0 = dfes->DofToVDof(local_node_num,0);
      id1 = dfes->DofToVDof(local_node_num,1);
      if (dim == 3) id2 = dfes->DofToVDof(local_node_num,2);

      x = dn_gf(id0);
      y = dn_gf(id1);
      if (dim == 3) z = dn_gf(id2);
   }  


   #ifdef MFEM_USE_MPI
   // communication to get to root proc
   MPI_Allreduce(&x,&x_glob,1,MPI_DOUBLE,MPI_SUM,dfes->GetComm());
   MPI_Allreduce(&y,&y_glob,1,MPI_DOUBLE,MPI_SUM,dfes->GetComm());
   if (dim == 3) MPI_Allreduce(&z,&z_glob,1,MPI_DOUBLE,MPI_SUM,dfes->GetComm());

   disp(0) = x_glob;
   disp(1) = y_glob;
   if (dim == 3) disp(2) = z_glob;
   #else
   disp(0) = x;
   disp(1) = y;
   if (dim == 3) disp(2) = z;
   #endif
}

void ImmersedSolid::PrintInfo()
{
   #ifdef MFEM_USE_MPI
   int fes_size0 = dfes->GlobalTrueVSize();
   int numelements = pmesh->GetNE();
   MPI_Allreduce(MPI_IN_PLACE, &numelements, 1, MPI_INT, MPI_SUM, vfes->GetComm());
   int global_elem = numelements;
   #else
   int fes_size0 = dfes->GetTrueVSize();
   int global_elem = pmesh->GetNE();
   #endif

   if (my_rank == 0)
   {
      mfem::out << "ImmersedSolid version: " << IMMERSEDSOLID_VERSION << std::endl
                  << "MFEM version: " << MFEM_VERSION << std::endl
                  << "Number of elements: " << global_elem << std::endl
                  << "Displacement #DOFs: " << fes_size0 << std::endl;
   }
}

#ifdef MFEM_USE_MPI
void ImmersedSolid::AddTractionBC(VectorCoefficient *coeff, Array<int> &attr)
{
   // note this is expected to be defined with respect to reference configuration
   traction_bcs.emplace_back(attr, coeff);

   if (verbose && pmesh->GetMyRank() == 0)
   {
      mfem::out << "Adding Immersed Traction BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }

   for (int i = 0; i < attr.Size(); ++i)
   {
      MFEM_ASSERT((traction_attr[i] && attr[i]) == 0,
                  "Duplicate boundary definition detected.");
      if (attr[i] == 1)
      {
         traction_attr[i] = 1;
      }
   }
}

void ImmersedSolid::AddTractionBC(VecFuncT *f, Array<int> &attr)
{
   // note this is expected to be defined with respect to reference configuration
   AddTractionBC(new VectorFunctionCoefficient(pmesh->Dimension(), f), attr);
}

void ImmersedSolid::AddFollowerPressureLoad(Coefficient *coeff, Array<int> &attr)
{
   // note this is expected to be defined with respect to current configuration
   // no minus sign is included. i.e t = p n. for a compressive load, specify p = -100 for example.
   
   MFEM_ASSERT(pmesh->Dimension() == 3,"Error: TO DO, follower pressure load only implemented in 3D");
   
   fp_bcs.emplace_back(attr, coeff);

   if (verbose && pmesh->GetMyRank() == 0)
   {
      mfem::out << "Adding Immersed Follower Pressure Load BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }

   for (int i = 0; i < attr.Size(); ++i)
   {
      MFEM_ASSERT((traction_attr[i] && attr[i]) == 0,
                  "Duplicate boundary definition detected.");
      if (attr[i] == 1)
      {
         traction_attr[i] = 1;
      }
   }
}

void ImmersedSolid::AddFollowerPressureLoad(ScalarFuncT *f, Array<int> &attr)
{
   AddFollowerPressureLoad(new FunctionCoefficient(f), attr);
}
#endif

void ImmersedSolid::writeRestarts(std::ofstream &sout, std::string solidmname, std::string dname, std::string ddotname, std::string d2dotname)
{
   // write solid mesh
   sout.open(solidmname.c_str(),std::ios::out);
   sout.precision(20);
   #ifdef MFEM_USE_MPI
   pmesh->ParPrint(sout);
   #else
   pmesh->Print(sout);
   #endif
   sout.close();

   // write solid state
   sout.open(dname.c_str(),std::ios::out);
   sout.precision(20);
   dn_gf.Save(sout);
   sout.close();

   sout.open(ddotname.c_str(),std::ios::out);
   sout.precision(20);
   ddotn_gf.Save(sout);
   sout.close();

   sout.open(d2dotname.c_str(),std::ios::out);
   sout.precision(20);
   d2dotn_gf.Save(sout);
   sout.close();
}

void ImmersedSolid::readRestarts(std::ifstream &in, BaseKDTreeNodalProjection* map_d, BaseKDTreeNodalProjection* map_ddot, BaseKDTreeNodalProjection* map_d2dot, std::string solidmname, std::string dname, std::string ddotname, std::string d2dotname)
{
   // read solid mesh
   Mesh lmesh;
   in.open(solidmname.c_str(),std::ios::in);
   lmesh.Load(in);
   in.close();

   in.open(dname.c_str(),std::ios::in);
   GridFunction d_in_gf(&lmesh,in);
   in.close();

   if (formulation_type == form_type::UpdatedLagrangian)
   {
      // pull back solid mesh to reference config
      // since we will be mapping to the initial mesh
      lmesh.EnsureNodes();
      GridFunction disp_gf(d_in_gf);
      disp_gf *= -1.0;
      lmesh.MoveNodes(disp_gf); // x0 = x - d
   }   

   in.open(ddotname.c_str(),std::ios::in);
   GridFunction ddot_in_gf(&lmesh,in);
   in.close();

   in.open(d2dotname.c_str(),std::ios::in);
   GridFunction d2dot_in_gf(&lmesh,in);
   in.close();

   double tol = 1e-10;
   map_d->Project(d_in_gf,tol);
   map_ddot->Project(ddot_in_gf,tol);
   map_d2dot->Project(d2dot_in_gf,tol);
}

void ImmersedSolid::UpdateSolidStateFromRestart()
{
   // update solid state after reading restart files
   #ifdef MFEM_USE_MPI
   dn_gf.ParallelProject(dn);
   ddotn_gf.ParallelProject(ddotn);
   d2dotn_gf.ParallelProject(d2dotn);
   #else
   dn_gf.GetTrueDofs(dn);
   ddotn_gf.GetTrueDofs(ddotn);
   d2dotn_gf.GetTrueDofs(d2dotn);
   #endif

   // update solid mesh (if UL formulation) //
   MoveMesh(dn_gf);
   UpdatePointSearch(dn_gf,dn);
   if (stabs_scaling)
   {
      #ifdef MFEM_USE_GSLIB
      finder_fluid->TransferIndicator(solid_bdr_nodes,*s_gf_fluid,dfes->GetOrdering());
      #endif
   }
}

void ImmersedSolid::BuildInitTransformation(ElementTransformation &T, const IntegrationPoint &ip)
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
         auto dof_trans = dfes->GetBdrElementVDofs(T.ElementNo,vdofs);
      }
      else if (T.ElementType == ElementTransformation::ELEMENT)
      {
         // element transformation
         auto dof_trans = dfes->GetElementVDofs(T.ElementNo,vdofs);
      }
      else
      {
         MFEM_ABORT("TO DO: ERROR, rebuilding init_transformation, not implemented element type");
      }
      
      orig_nodes.GetSubVector(vdofs,vals);
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

GridFunction *ImmersedSolid::ComputeNodallyAveragedStress(nodal_avg::StressQuantity measure, GridFunction &d_gf)
{
   const FiniteElement *fe;
   ElementTransformation *Tr;
   DenseMatrix gradu;
   int intorder;
   int dim = d_gf.FESpace()->GetMesh()->SpaceDimension();
   Vector u(dim);
   u = 0.0;
   gradu.SetSize(dim);
   gradu = 0.0;

   DenseMatrix stress(dim);
   stress = 0.0;

   Vector shape, local_dat;
   DenseMatrix dshape, dshapedx, EF;

   Array<int> vdofs, vdofs_stress;

   if (nodal_stress_fec == nullptr)
   {
      // create a linear fes for stress output
      nodal_stress_fec = new H1_FECollection(orderd, dim);
      #ifdef MFEM_USE_MPI
      nodal_stress_fes = new ParFiniteElementSpace(pmesh, nodal_stress_fec);
      #else
      nodal_stress_fes = new FiniteElementSpace(pmesh, nodal_stress_fec);
      #endif
   }

   // get or initialize the target field for this measure
   NodalField &nodal_field = nodal_stress_fields[measure];
   if (nodal_field.FESpace() == nullptr)
   {
      nodal_field.SetSpace(nodal_stress_fes);
   }
   nodal_field = 0.0;

   // AccumulateAndCountZones.
   Array<int> zones_per_vdof;
   zones_per_vdof.SetSize(nodal_field.FESpace()->GetVSize());
   zones_per_vdof = 0;

   for (int i = 0; i < d_gf.FESpace()->GetNE(); i++)
   {
      #ifdef MFEM_USE_MPI
      auto par_d_gf = dynamic_cast<ParGridFunction *>(&d_gf);
      fe = par_d_gf->ParFESpace()->GetFE(i);
      #else
      fe = d_gf.FESpace()->GetFE(i);
      #endif
      Tr = d_gf.FESpace()->GetElementTransformation(i);
      d_gf.FESpace()->GetElementVDofs(i,vdofs);
      nodal_field.FESpace()->GetElementVDofs(i,vdofs_stress);

      int nd = fe->GetDof();

      shape.SetSize(nd);
      dshape.SetSize(nd, dim);
      dshapedx.SetSize(nd, dim);

      local_dat.SetSize(nd*dim);
      d_gf.GetSubVector(vdofs,local_dat);
      EF.UseExternalData(local_dat.GetData(), nd, dim);

      // loop over nodes
      const IntegrationRule *ir = &fe->GetNodes();
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir->IntPoint(j);
         Tr->SetIntPoint(&ip);

         fe->CalcShape(ip, shape);
         fe->CalcDShape(ip, dshape);

         // compute dshapedx using appropriate transformation for computing du/dX and F = dx/dX
         if (formulation_type == form_type::TotalLagrangian)
         {
            Mult(dshape, Tr->InverseJacobian(), dshapedx); // (dN/dX)
         }
         else if (formulation_type == form_type::UpdatedLagrangian)
         {
            BuildInitTransformation(*Tr,ip);
            Mult(dshape, init_trans->InverseJacobian(), dshapedx); // (dN/dX)
         }
         else
         {
            MFEM_ABORT("Error: Unsupported formulation type");
         }

         MultAtB(EF, dshapedx, gradu); // dudX
         EF.MultTranspose(shape, u); // u
         
         // compute strain measures, material properties at the current integration point, stored for all operations later on
         if (formulation_type == form_type::TotalLagrangian)
         {
            mat_type->SetStrainMeasures(u,gradu,ip,*Tr);
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
         mat_type->EvalCauchyStress(stress);

         double stress_val = StressMeasureValue(stress, measure, dim);

         nodal_field(vdofs_stress[j]) += stress_val;
         zones_per_vdof[vdofs_stress[j]]++;
      }
   }

   #ifdef MFEM_USE_MPI
   // average across processors
   // Count the zones globally.
   GroupCommunicator &gcomm = nodal_field.ParFESpace()->GroupComm();
   gcomm.Reduce<int>(zones_per_vdof, GroupCommunicator::Sum);
   gcomm.Bcast(zones_per_vdof);

   // Accumulate for all vdofs.
   gcomm.Reduce<double>(nodal_field.GetData(), GroupCommunicator::Sum);
   gcomm.Bcast<double>(nodal_field.GetData());
   #endif

   // Compute means.
   for (int i = 0; i < nodal_field.Size(); i++)
   {
      const int nz = zones_per_vdof[i];
      if (nz)
      {
         nodal_field(i) /= nz;
      }
   }

   return &nodal_field;
}

GridFunction *ImmersedSolid::ComputeNodallyAveragedStrain(nodal_avg::StrainQuantity measure, GridFunction &d_gf)
{
   const FiniteElement *fe;
   ElementTransformation *Tr;
   DenseMatrix gradu;
   int intorder;
   int dim = d_gf.FESpace()->GetMesh()->SpaceDimension();
   Vector u(dim);
   u = 0.0;
   gradu.SetSize(dim);
   gradu = 0.0;

   DenseMatrix strain(dim);
   strain = 0.0;

   Vector shape, local_dat;
   DenseMatrix dshape, dshapedx, EF;

   Array<int> vdofs, vdofs_strain;

   if (nodal_strain_fec == nullptr)
   {
      // create a linear fes for strain output
      nodal_strain_fec = new H1_FECollection(orderd, dim);
      #ifdef MFEM_USE_MPI
      nodal_strain_fes = new ParFiniteElementSpace(pmesh, nodal_strain_fec);
      #else
      nodal_strain_fes = new FiniteElementSpace(pmesh, nodal_strain_fec);
      #endif
   }

   // get or initialize the target field for this measure
   NodalField &nodal_field = nodal_strain_fields[measure];
   if (nodal_field.FESpace() == nullptr)
   {
      nodal_field.SetSpace(nodal_strain_fes);
   }
   nodal_field = 0.0;

   // AccumulateAndCountZones.
   Array<int> zones_per_vdof;
   zones_per_vdof.SetSize(nodal_field.FESpace()->GetVSize());
   zones_per_vdof = 0;

   for (int i = 0; i < d_gf.FESpace()->GetNE(); i++)
   {
      #ifdef MFEM_USE_MPI
      auto par_d_gf = dynamic_cast<ParGridFunction *>(&d_gf);
      fe = par_d_gf->ParFESpace()->GetFE(i);
      #else
      fe = d_gf.FESpace()->GetFE(i);
      #endif
      Tr = d_gf.FESpace()->GetElementTransformation(i);
      d_gf.FESpace()->GetElementVDofs(i,vdofs);
      nodal_field.FESpace()->GetElementVDofs(i,vdofs_strain);

      int nd = fe->GetDof();

      shape.SetSize(nd);
      dshape.SetSize(nd, dim);
      dshapedx.SetSize(nd, dim);

      local_dat.SetSize(nd*dim);
      d_gf.GetSubVector(vdofs,local_dat);
      EF.UseExternalData(local_dat.GetData(), nd, dim);

      // loop over nodes
      const IntegrationRule *ir = &fe->GetNodes();
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir->IntPoint(j);
         Tr->SetIntPoint(&ip);

         fe->CalcShape(ip, shape);
         fe->CalcDShape(ip, dshape);

         // compute dshapedx using appropriate transformation for computing du/dX and F = dx/dX
         if (formulation_type == form_type::TotalLagrangian)
         {
            Mult(dshape, Tr->InverseJacobian(), dshapedx); // (dN/dX)
         }
         else if (formulation_type == form_type::UpdatedLagrangian)
         {
            BuildInitTransformation(*Tr,ip);
            Mult(dshape, init_trans->InverseJacobian(), dshapedx); // (dN/dX)
         }
         else
         {
            MFEM_ABORT("Error: Unsupported formulation type");
         }

         MultAtB(EF, dshapedx, gradu); // dudX
         EF.MultTranspose(shape, u); // u
         
         // compute strain measures, material properties at the current integration point, stored for all operations later on
         if (formulation_type == form_type::TotalLagrangian)
         {
            mat_type->SetStrainMeasures(u,gradu,ip,*Tr);
         }
         else if (formulation_type == form_type::UpdatedLagrangian)
         {
            mat_type->SetStrainMeasures(u,gradu,ip,*init_trans);
         }
         else
         {
            MFEM_ABORT("Error: Unsupported formulation type");
         }      

         // compute strain
         strain = mat_type->GetE();

         double strain_val = StrainMeasureValue(strain, measure, dim);

         nodal_field(vdofs_strain[j]) += strain_val;
         zones_per_vdof[vdofs_strain[j]]++;
      }
   }

   #ifdef MFEM_USE_MPI
   // average across processors
   // Count the zones globally.
   GroupCommunicator &gcomm = nodal_field.ParFESpace()->GroupComm();
   gcomm.Reduce<int>(zones_per_vdof, GroupCommunicator::Sum);
   gcomm.Bcast(zones_per_vdof);

   // Accumulate for all vdofs.
   gcomm.Reduce<double>(nodal_field.GetData(), GroupCommunicator::Sum);
   gcomm.Bcast<double>(nodal_field.GetData());
   #endif

   // Compute means.
   for (int i = 0; i < nodal_field.Size(); i++)
   {
      const int nz = zones_per_vdof[i];
      if (nz)
      {
         nodal_field(i) /= nz;
      }
   }

   return &nodal_field;
}

ImmersedSolid::~ImmersedSolid()
{
   #ifdef MFEM_USE_GSLIB
   if (finder_fluid) finder_fluid->FreeData();
   delete finder_fluid;
   #endif
   
   delete dfec;
   delete dfes;

   delete density_diff;
   delete iSolidNLF;

   delete KCt;
   delete CKCt;

   #ifdef MFEM_USE_MPI
   delete Ft_form;
   #endif

   delete init_trans;
   delete nodal_stress_fec;
   delete nodal_stress_fes;
   delete nodal_strain_fec;
   delete nodal_strain_fes;
   
}