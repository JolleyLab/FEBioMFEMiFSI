#include "MFEMImmersedFEBioSolid.h"

#include <FECore/FEDomain.h>
#include <FECore/FEModelLoad.h>
#include <FEBioMech/FEElasticDomain.h>

#include "../MFEMImSolid/MFEMImmersedElasticDomain.h"
#include "mfem.hpp"
#include <general/forall.hpp>



using namespace mfem;
using namespace immersed_solid;

#ifdef MFEM_USE_MPI
ImmersedFEBioSolid::ImmersedFEBioSolid(FEModel* fem_, mfem::ParMesh *solid_mesh, int order_, mfem::Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
   mfem::ParFiniteElementSpace *vfes_, mfem::ParMesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_)
: ImmersedSolid(solid_mesh, order_, init_rho_s_, fluid_density_, mu_f_, accel_, vfes_, fluid_mesh, mat_, orderingtype_, form_type_, subtract_visc_), fem(fem_)
#else
ImmersedFEBioSolid::ImmersedFEBioSolid(FEModel* fem_, mfem::Mesh *solid_mesh, int order_, mfem::Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
   mfem::FiniteElementSpace *vfes_, mfem::Mesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_)
: ImmersedSolid(solid_mesh, order_, init_rho_s_, fluid_density_, mu_f_, accel_, vfes_, fluid_mesh, mat_, orderingtype_, form_type_, subtract_visc_), fem(fem_)
#endif
{
   int dfes_truevsize = dfes->GetTrueVSize();
   int dfes_vsize = dfes->GetVSize();

   // setup and initialize additional variables
   dnp1_gf.SetSpace(dfes);
   dnp1_gf = 0.0;
   ddotnp1_gf.SetSpace(dfes);
   ddotnp1_gf = 0.0;
   d2dotnp1_gf.SetSpace(dfes);
   d2dotnp1_gf = 0.0;
   
   // setup residual vector
   R.resize(dfes_vsize,0.0);
   m_R = new MFEMGlobalVector(fem,dfes_vsize,R,Fr);

   // setup linear system to assemble FEBio into
   CRSSparseMatrix *pS = new CRSSparseMatrix(0);
   if (m_K) delete m_K;
   m_K = new MFEMGlobalMatrix(pS); // ps owned by MFEMGlobalMatrix

   #ifdef MFEM_USE_MPI
   pK.Clear();
   pK.SetType(mfem::Operator::Hypre_ParCSR);
   #endif
}


void ImmersedFEBioSolid::Setup(double dt, double alpha_m_, double alpha_f_, double gamma_, double beta_)
{
   if (my_rank == 0)
   {
      mfem::out << "Setup Immersed Solid" << std::endl;     
   }
   PrintInfo();

   this->dt = dt;
   
   // // setup point search
   // #ifdef MFEM_USE_GSLIB   

   // #ifdef MFEM_USE_MPI
   // finder_fluid = new ImmersedFSI_FindPts(vfes->GetComm());
   // #else
   // finder_fluid = new ImmersedFSI_FindPts();
   // #endif

   // finder_fluid->Setup(*fluid_pmesh);

   // #ifdef MFEM_USE_MPI
   // finder_fluid->SetupLocalDofToMasterRankMap(vfes);
   // finder_fluid->SetupLocalTDofToLDofSolidMap(dfes);
   // #endif
   // #endif

   // true_interp_points_in_solid.SetSize(dfes->GetTrueVSize());
   // true_interp_points_in_solid = 0.0;

   // interp_points_in_solid = *pmesh->GetNodes();
   // mfem::GridFunction *nodes = pmesh->GetNodes();
   // #ifdef MFEM_USE_MPI
   // dynamic_cast<mfem::ParGridFunction *>(nodes)->ParallelProject(true_interp_points_in_solid);
   // #else
   // nodes->GetTrueDofs(true_interp_points_in_solid);
   // #endif

   // CheckMeshIntersection(true_interp_points_in_solid);

   beta = beta_;
   gamma = gamma_;
   alpha_m = alpha_m_;
   alpha_f = alpha_f_;

   // set equivalent FEBio data
   FETimeInfo& tp = fem->GetTime();
   tp.timeIncrement = dt;
   tp.alpham = alpha_m;
   tp.alphaf = alpha_f;
   tp.beta = beta;
   tp.gamma = gamma;

   // setup vector of body forces
   for (int i = 0; i < fem->ModelLoads(); ++i)
	{
      // body loads
      MFEMBodyForce* bl = dynamic_cast<MFEMBodyForce*>(fem->ModelLoad(i));
      if (bl && bl->IsActive())
      {  
         m_body_forces.push_back(bl);        
      }
	}

   // setup FEBio domain
   // FEMesh &mesh = fem->GetMesh();
   // for (int nd = 0; nd<mesh.Domains(); ++nd)
   // {
   //    FEDomain& dom = mesh.Domain(nd);
   //    MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
   //    if (mfd)
   //    {
   //       mfd->SetFESpace(dfes);
   //       mfd->SetFluidProperties(fluid_density,mu_f);
   //       #ifdef MFEM_USE_MPI
   //       mfd->SetAdditionalData(subtract_fluid_visc);
   //       #else
   //       mfd->SetAdditionalData(subtract_fluid_visc,dynamic_cast<mfem::SubMesh*>(pmesh),ddot_global_gf->FESpace());
   //       #endif
   //    }			
   // }
   SetupFEDomains();

   // set initial conditions
   #ifdef MFEM_USE_MPI
   dn_gf.ParallelProject(dn);
   #else
   dn_gf.GetTrueDofs(dn);
   #endif

   // // setup VMS stab scaling
   // if (stabs_scaling)
   // {
   //    // setup interpolation points (only boundary dofs of the solid)
   //    mfem::Array<int> bdr_attr(pmesh->bdr_attributes.Max());
   //    bdr_attr = 1;
   //    dfes->GetEssentialTrueDofs(bdr_attr,solid_bdr_tdofs);
   //    solid_bdr_nodes.SetSize(solid_bdr_tdofs.Size());

   //    const int Ndofs = solid_bdr_tdofs.Size();
   //    const auto id = solid_bdr_tdofs.Read();
   //    auto FV_nodes = true_interp_points_in_solid.Read();
   //    auto FV = solid_bdr_nodes.ReadWrite();
   //    mfem::forall(Ndofs, [=] MFEM_HOST_DEVICE (int i) { FV[i] = FV_nodes[id[i]]; });
   // }
}

void ImmersedFEBioSolid::SetupFEDomains()
{
   // setup FEBio domain
   FEMesh &mesh = fem->GetMesh();
   for (int nd = 0; nd<mesh.Domains(); ++nd)
   {
      FEDomain& dom = mesh.Domain(nd);
      MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (mfd)
      {
         mfd->SetFESpace(dfes);
         mfd->SetFluidProperties(fluid_density,mu_f);
         #ifdef MFEM_USE_MPI
         mfd->SetAdditionalData(subtract_fluid_visc);
         #else
         mfd->SetAdditionalData(subtract_fluid_visc,dynamic_cast<mfem::SubMesh*>(pmesh),ddot_global_gf->FESpace());
         #endif
      }			
   }
}

void ImmersedFEBioSolid::SetSubMeshMaps(mfem::TransferMap *submeshsolid_to_global_d_, mfem::TransferMap* global_to_submeshsolid_d_, mfem::GridFunction *ddot_global_gf_)
{
   submeshsolid_to_global_d = submeshsolid_to_global_d_;
   global_to_submeshsolid_d = global_to_submeshsolid_d_;
   ddot_global_gf = ddot_global_gf_;
   d2dot_global_gf = new mfem::GridFunction(*ddot_global_gf);
   *d2dot_global_gf = 0.0;
   dnp1_global_gf = new mfem::GridFunction(*ddot_global_gf);
   *dnp1_global_gf = 0.0;

   delvdot_global_gf.SetSpace(ddot_global_gf->FESpace());
   delvdot_global_gf = 0.0;
}

void ImmersedFEBioSolid::Predictors(mfem::Vector &vf, mfem::Vector &vdotf)
{
   // set ddotnp1 and d2dotnp1 from fluid predictors for consistency!
   // ddotnp1 = v_np1
   temp_del_gf.SetFromTrueDofs(vf);
   #ifdef MFEM_USE_GSLIB
   finder_fluid->Interpolate(interp_points_in_solid,temp_del_gf,ddotnp1_gf,orderingtype);
   #endif
   #ifdef MFEM_USE_MPI
   ddotnp1_gf.ParallelProject(ddotnp1);
   #else
   ddotnp1_gf.GetTrueDofs(ddotnp1);
   #endif

   // d2dotnp1 = vdot_np1
   temp_del_gf.SetFromTrueDofs(vdotf);
   #ifdef MFEM_USE_GSLIB
   finder_fluid->Interpolate(interp_points_in_solid,temp_del_gf,d2dotnp1_gf,orderingtype);
   #endif
   #ifdef MFEM_USE_MPI
   d2dotnp1_gf.ParallelProject(d2dotnp1);
   #else
   d2dotnp1_gf.GetTrueDofs(d2dotnp1);   
   #endif

   #ifndef MFEM_USE_MPI
   // transfer gf defined on solid submesh to global mesh (FEBio representation)
   *ddot_global_gf = 0.0; // not sure if needed?
   *d2dot_global_gf = 0.0;
   submeshsolid_to_global_d->Transfer(ddotnp1_gf,*ddot_global_gf);
   submeshsolid_to_global_d->Transfer(d2dotnp1_gf,*d2dot_global_gf);
   mfem::Vector &ddot = *ddot_global_gf;
   mfem::Vector &d2dot = *d2dot_global_gf;
   mfem::Vector &d = *dnp1_global_gf;
   #else
   mfem::Vector &ddot = ddotnp1_gf;
   mfem::Vector &d2dot = d2dotnp1_gf;
   mfem::Vector &d = dnp1_gf;
   #endif
   int dofVIMX = fem->GetDOFIndex("im_vx");
   int dofVIMY = fem->GetDOFIndex("im_vy");
   int dofVIMZ = fem->GetDOFIndex("im_vz");   

   // Update nodal data, which stores values at n and np1
   FETimeInfo& tp = fem->GetTime();
	FEMesh& mesh = fem->GetMesh();
   int NN = mesh.Nodes();
	for (int i=0; i<NN; ++i)
	{      
		FENode& ni = mesh.Node(i);
      // store prev state
		ni.m_rp = ni.m_rt;
		ni.m_vp = ni.get_vec3d(dofVIMX, dofVIMY, dofVIMZ);
		ni.m_ap = ni.m_at;
		ni.UpdateValues();

      // solid predictors
      int idx = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 0*NN : i*3 + 0;
      int idy = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 1*NN : i*3 + 1;
      int idz = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 2*NN : i*3 + 2;

      // acceleration
      ni.m_at(0) = d2dot[idx];
      ni.m_at(1) = d2dot[idy];
      ni.m_at(2) = d2dot[idz];

      // velocity
      vec3d vs(ddot[idx],ddot[idy],ddot[idz]);
      ni.set_vec3d(dofVIMX, dofVIMY, dofVIMZ, vs);

      // position from displacement predictor
      double fac_n = 0.5*dt*dt*(1.0 - 2.0*beta);
      double fac_np1 = 0.5*dt*dt*(2.0*beta);
      vec3d dn(ni.m_rt(0), ni.m_rt(1), ni.m_rt(2));
      dn -= ni.m_r0; // compute dn = xt - x0
      double dnp1_x = dn(0) + dt*ni.m_vp(0) + fac_n*ni.m_ap(0) + fac_np1*ni.m_at(0);
      double dnp1_y = dn(1) + dt*ni.m_vp(1) + fac_n*ni.m_ap(1) + fac_np1*ni.m_at(1);
      double dnp1_z = dn(2) + dt*ni.m_vp(2) + fac_n*ni.m_ap(2) + fac_np1*ni.m_at(2);
      ni.m_rt(0) = dnp1_x + ni.m_r0(0);
      ni.m_rt(1) = dnp1_y + ni.m_r0(1);
      ni.m_rt(2) = dnp1_z + ni.m_r0(2);

      // set mfem dnp1
      d[idx] = dnp1_x;
      d[idy] = dnp1_y;
      d[idz] = dnp1_z;
   }

   #ifndef MFEM_USE_MPI
   // Transfer global dnp1 to submesh
   dnp1_gf = 0.0; // not sure if needed?
   global_to_submeshsolid_d->Transfer(*dnp1_global_gf,dnp1_gf);
   dnp1_gf.GetTrueDofs(dnp1);
   #else
   dnp1_gf.ParallelProject(dnp1);
   #endif

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
   for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
      MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (mfd)
      {
         mfd->PreSolveUpdate(tp);
      }	
	}

   // store tn position of solid
   interp_points_in_solid_tn = interp_points_in_solid;
   if (stabs_scaling) solid_bdr_nodes_tn = solid_bdr_nodes;
}

void ImmersedFEBioSolid::Initiator()
{
   // update element data (material point data, which is data stored
   // at the integration points) to the generalized midpoint
	fem->Update();

   // update surface data
   if (immersed_surface_loads)
   {
      FEMesh &mesh = fem->GetMesh();
      FETimeInfo &tp = fem->GetTime();
      for (int i=0; i<mesh.Surfaces(); ++i)
      {
         FESurface* ps = &mesh.Surface(i);
         // ps->Update(tp);
         for (int e=0; e<ps->Elements(); ++e) 
         {
            FESurfaceElement& el = ps->Element(e);

            // check if active
            if (el.GetID() == -1) continue;         

            int nint = el.GaussPoints();
            int neln = el.Nodes();

            vec3d rt[FEElement::MAX_NODES];
            // ps->NodalCoordinates(el, rt);
            for (int ii = 0; ii < neln; ++ii) rt[ii] = mesh.Node(el.m_lnode[ii]).m_rt;  

            // // ps->GetNodalCoordinates(el, alpha_f, rt);
            // for (int j = 0; j<neln; ++j) 
            // {
            //    FENode& node = mesh.Node(el.m_lnode[j]);
            //    rt[j] = node.m_rt*alpha_f + node.m_rp*(1.0 - alpha_f);
            // }
                  

            for (int n = 0; n < nint; ++n)
            {
               FESurfaceMaterialPoint& mp = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));

               // material point coordinates
               mp.m_rt = el.Evaluate(rt, n); // sets current position at int point to midpoint value

               double* Gr = el.Gr(n);
               double* Gs = el.Gs(n);

               mp.dxr = vec3d(0, 0, 0);
               mp.dxs = vec3d(0, 0, 0);
               for (int a = 0; a < neln; ++a)
               {
                  mp.dxr += rt[a] * Gr[a];
                  mp.dxs += rt[a] * Gs[a];
               }
            }
         }
      }
   }

	for (int i = 0; i < fem->ModelLoads(); ++i)
	{
		FEModelLoad* pl = fem->ModelLoad(i);
		if (pl->IsActive()) pl->PrepStep();
	}

   // need to update dg so that I can move the mfem mesh and transfer indicator below
   // move mesh and reset point search
   add((1.0 - alpha_f),dn,alpha_f,dnp1,dg);
   dtemp_gf.SetFromTrueDofs(dg);
   MoveMesh(dtemp_gf); // contains current displacements
   // UpdatePointSearch(dtemp_gf,dg); // not needed since motion map for lambda terms fixed at tn

   if (stabs_scaling)
   {
      #ifdef MFEM_USE_GSLIB
      MFEM_ASSERT(dfes->GetOrdering() == s_gf_fluid->FESpace()->GetOrdering(),"Ordering of solid and fluid gridfunctions must match for special VMS stabs scaling");
      
      // finder_fluid->TransferIndicator(solid_bdr_nodes,*s_gf_fluid,dfes->GetOrdering());
      finder_fluid->TransferIndicator(solid_bdr_nodes_tn,*s_gf_fluid,dfes->GetOrdering());
      // note internally s_gf_fluid is zeroed prior to interpolation of s_gf to its fluid counterpart
      #endif
   }
}

void ImmersedFEBioSolid::Correctors(mfem::Vector &del)
{
   // project delvdot to solid mesh and update data
   temp_del_gf.SetFromTrueDofs(del);
   #ifdef MFEM_USE_GSLIB
   // finder_fluid->Interpolate(interp_points_in_solid,temp_del_gf,temp_delvdot_ib_gf,orderingtype);
   finder_fluid->Interpolate(interp_points_in_solid_tn,temp_del_gf,temp_delvdot_ib_gf,orderingtype);
   #endif
   #ifdef MFEM_USE_MPI
   temp_delvdot_ib_gf.ParallelProject(delvdot_ib);
   #else
   temp_delvdot_ib_gf.GetTrueDofs(delvdot_ib);
   #endif

   d2dotnp1.Add(-1.0,delvdot_ib);
   ddotnp1.Add(-gamma*dt,delvdot_ib); // only updated to keep ddotn_gf current
   dnp1.Add(-beta*dt*dt,delvdot_ib); // for moving mfem mesh

   ////////////////////////////////////////////////////////////////

   #ifndef MFEM_USE_MPI
   // transfer submesh delvdot to global
   delvdot_global_gf = 0.0;
   submeshsolid_to_global_d->Transfer(temp_delvdot_ib_gf,delvdot_global_gf);
   mfem::Vector &delvdot = delvdot_global_gf;

   // transfer displacement to global
   dnp1_gf.SetFromTrueDofs(dnp1);
   *dnp1_global_gf = 0.0;
   submeshsolid_to_global_d->Transfer(dnp1_gf,*dnp1_global_gf);
   mfem::Vector &d = *dnp1_global_gf;
   #else
   mfem::Vector &delvdot = temp_delvdot_ib_gf;
   dnp1_gf.SetFromTrueDofs(dnp1);
   mfem::Vector &d = dnp1_gf;
   #endif
   
   // update nodal data
	FEMesh& mesh = fem->GetMesh();
   int dofVIMX = fem->GetDOFIndex("im_vx");
   int dofVIMY = fem->GetDOFIndex("im_vy");
   int dofVIMZ = fem->GetDOFIndex("im_vz");
   int NN = mesh.Nodes();
	for (int i=0; i<NN; ++i)
	{      
		FENode& ni = mesh.Node(i);
      
      int idx = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 0*NN : i*3 + 0;
      int idy = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 1*NN : i*3 + 1;
      int idz = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 2*NN : i*3 + 2;

      // acceleration
      ni.m_at(0) -= delvdot[idx];
      ni.m_at(1) -= delvdot[idy];
      ni.m_at(2) -= delvdot[idz];

      // velocity
      vec3d vt = ni.get_vec3d(dofVIMX, dofVIMY, dofVIMZ);
      vt(0) -= gamma*dt*delvdot[idx];
      vt(1) -= gamma*dt*delvdot[idy];
      vt(2) -= gamma*dt*delvdot[idz];
      ni.set_vec3d(dofVIMX, dofVIMY, dofVIMZ, vt);

      // position from updated displacement
      ni.m_rt(0) = d[idx] + ni.m_r0(0);
      ni.m_rt(1) = d[idy] + ni.m_r0(1);
      ni.m_rt(2) = d[idz] + ni.m_r0(2);
   }
}

void ImmersedFEBioSolid::Update()
{
   // update dn_gf and ddotn_gf in mfem,
   dn = dnp1;
   dn_gf.SetFromTrueDofs(dn);

   // update ddotn_gf in mfem
   ddotn = ddotnp1;
   ddotn_gf.SetFromTrueDofs(ddotn);

   // update d2dotn_gf in mfem
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

void ImmersedFEBioSolid::ComputeSolidResidualAndSpread2Fluid(mfem::Vector &lambda_fsi_vec)
{
   sw_res->Start();

   // Assemble residual vector from FEBio Solid Domain
   // get the time information
	const FETimeInfo& tp = fem->GetTime();

	// setup the global vector
	zero(R);

	// calculate residual

   // domain residual and body forces
   auto copy_fes = dfes;
   FEMesh& mesh = fem->GetMesh();
   for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
      MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (mfd)
      {
         mfd->AssembleDomainResidual(*m_R,m_body_forces);
      }	
	}

   // surface forces
   if (immersed_surface_loads)
   {
      auto solid_submesh = dynamic_cast<mfem::SubMesh*>(pmesh);
      FESurface& surf = mesh.Surface(0); // assumes one surface was built for all loads
      int NBE = surf.Elements();
      int dofPerNode = dfes->GetVDim();
      mfem::Array<int> vdofs, vdofs_febio;

      double dt = tp.timeIncrement;
      double am = tp.alpham;
      double av = tp.alphaf*tp.gamma*dt;
      double af = tp.alphaf*tp.beta*dt*dt;
      std::vector<double> factors({am, av, af, dt});

      for (int b=0; b<NBE; ++b)
      {
         const int bdr_attr = pmesh->GetBdrAttribute(b);
         for (int i=0; i < fem->ModelLoads(); ++i)
         {
            FEImmersedSurfaceLoad* sl = dynamic_cast<FEImmersedSurfaceLoad*>(fem->ModelLoad(i));
            if (sl && sl->IsActive())
            {
               int curr_bdr_attr = sl->GetBoundaryID();
               if (curr_bdr_attr != bdr_attr-1) { continue; }

               vector<double> fe;
               vec3d re[FEElement::MAX_NODES];
               std::vector<double> G(dofPerNode, 0.0);

               // get the next element
               FESurfaceElement& el = surf.Element(b);

               // get corresponding mfem element id
               int submesh_id = -1;
               if (solid_submesh)
               {
                  // serial
                  MFEM_ABORT("Not implemented");
               }         
               else
               {
                  // parallel
                  submesh_id = el.GetID() - 1;
               }

               // init the element vector
               int neln = el.ShapeFunctions(orderd);
               int ndof = dofPerNode * neln;
               fe.assign(ndof, 0.0);

               // get the nodal coordinates
               if (sl->IsLinear())
                  surf.GetReferenceNodalCoordinates(el, re);
               else
                  surf.GetNodalCoordinates(el, re);

               // calculate element vector
               FESurfaceDofShape dof_a;
               double* w = el.GaussWeights();
               int nint = el.GaussPoints();
               for (int n = 0; n < nint; ++n)
               {
                  FESurfaceMaterialPoint& pt = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));

                  // kinematics at integration points
                  pt.dxr = el.eval_deriv1(re, n);
                  pt.dxs = el.eval_deriv2(re, n);

                  pt.m_shape = el.H(n);

                  double* H = el.H(orderd, n);
                  double* Hr = el.Gr(orderd, n);
                  double* Hs = el.Gr(orderd, n);

                  // put it all together
                  for (int j = 0; j<neln; ++j)
                  {
                     // shape function and derivatives
                     dof_a.index = j;
                     dof_a.shape = H[j];
                     dof_a.shape_deriv_r = Hr[j];
                     dof_a.shape_deriv_s = Hs[j];

                     // evaluate the integrand
                     sl->ResidualIntegrand(pt, dof_a, G, factors);

                     for (int k = 0; k < dofPerNode; ++k)
                     {
                        fe[dofPerNode * j + k] -= G[k] * w[n];
                     }
                  }
               }

               // assemble element 'fe'-vector into global R vector
               dfes->GetBdrElementVDofs(submesh_id,vdofs);
               
               // reorder to match febio elem vector [xyz,xyz,xyz,...]
               vdofs_febio.SetSize(vdofs.Size());
               int nd = el.Nodes();
               for (int j=0; j < nd; ++j)
               {
                  vdofs_febio[3*j]   = vdofs[j];
                  vdofs_febio[3*j+1] = vdofs[j+nd];
                  vdofs_febio[3*j+2] = vdofs[j+2*nd];
               }

               m_R->Assemble(vdofs_febio,fe);
            }
         }
      }
   }

   // TODO: contact

   ////////////////////////////////////////////////////////
   // Get pointer to mfem::Vector version of m_R
   mfem::Vector &residual_vector = m_R->GetMfemVector();

   #ifdef MFEM_USE_MPI
   // Prolongate the residual vector
   auto P = dfes->GetProlongationMatrix();
   P->MultTranspose(residual_vector,res);
   mfem::Vector &res_vec = res;
   #else
   mfem::Vector &res_vec = residual_vector;
   #endif

   // if (disp_ess_tdof_list.Size() > 0)
   // {
   //    // zero essential displacement dofs
   //    const int N = disp_ess_tdof_list.Size();
   //    const auto idx = disp_ess_tdof_list.Read();
   //    auto r_RW = res_vec.ReadWrite();
   //    MFEM_FORALL(i, N, r_RW[idx[i]] = 0.0;);
   // }

   sw_res->Stop();

   sw_fsires->Start();

   // Spread to fluid nodes
   #ifdef MFEM_USE_GSLIB
   #ifdef MFEM_USE_MPI
   mfem::HypreParMatrix *C = dynamic_cast<mfem::HypreParMatrix *>(&finder_fluid->GetCouplingMatrix());
   #else
   mfem::SparseMatrix *C = dynamic_cast<mfem::SparseMatrix *>(&finder_fluid->GetCouplingMatrix());
   #endif
   #else
   #ifdef MFEM_USE_MPI
   mfem::HypreParMatrix *C = nullptr;
   #else
   mfem::SparseMatrix *C = nullptr;
   #endif
   #endif

   C->Mult(res_vec,lambda_fsi_vec);

   sw_fsires->Stop();
}

mfem::Operator &ImmersedFEBioSolid::FormLHS()
{
   sw_tan->Start();

   const FETimeInfo& tp = fem->GetTime();

	// get the mesh
	FEMesh& mesh = fem->GetMesh();

   // setup global matrix
   ReformStiffness();

   // form the stiffness matrix
   double dt = tp.timeIncrement;
   double am = tp.alpham;
   double av = tp.alphaf*tp.gamma*dt;
   double af = tp.alphaf*tp.beta*dt*dt;
   std::vector<double> factors({am, av, af, dt});
   for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
      MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (mfd)
      {
         mfd->AssembleDomainStiffness(*m_K,am,av,af,m_body_forces);
      }	
	}

   // surface forces
   if (immersed_surface_loads)
   {
      auto solid_submesh = dynamic_cast<mfem::SubMesh*>(pmesh);
      FESurface& surf = mesh.Surface(0); // assumes one surface was built for all loads
      int NBE = surf.Elements();

      FEElementMatrix ke;
      int dofPerNode_a = dfes->GetVDim();
      int dofPerNode_b = dfes->GetVDim();
      int order_a = orderd;
      int order_b = orderd;

      vec3d rt[FEElement::MAX_NODES];
      matrix kab(dofPerNode_a, dofPerNode_b);
      FESurfaceDofShape dof_a, dof_b;

      mfem::Array<int> vdofs, vdofs_febio;
      for (int b=0; b<NBE; ++b)
      {
         const int bdr_attr = pmesh->GetBdrAttribute(b);
         for (int i=0; i < fem->ModelLoads(); ++i)
         {
            FEImmersedSurfaceLoad* sl = dynamic_cast<FEImmersedSurfaceLoad*>(fem->ModelLoad(i));
            if (sl && sl->IsActive())
            {
               int curr_bdr_attr = sl->GetBoundaryID();
               // if (curr_bdr_attr != bdr_attr-1) { continue; }
               if (curr_bdr_attr-1 != bdr_attr-1) { continue; } // RB: 9-9-25

               // get the next element
               FESurfaceElement& el = surf.Element(b);

               // get corresponding mfem element id
               int submesh_id = -1;
               if (solid_submesh)
               {
                  // serial
                  MFEM_ABORT("Not implemented");
               }         
               else
               {
                  // parallel
                  submesh_id = el.GetID() - 1;
               }

               // get vdofs
               dfes->GetBdrElementVDofs(submesh_id,vdofs);
               // reorder to match febio elem vector [xyz,xyz,xyz,...]
               vdofs_febio.SetSize(vdofs.Size());
               int nd = el.Nodes();
               for (int j=0; j < nd; ++j)
               {
                  vdofs_febio[3*j]   = vdofs[j];
                  vdofs_febio[3*j+1] = vdofs[j+nd];
                  vdofs_febio[3*j+2] = vdofs[j+2*nd];
               }

               // shape functions
               int neln = el.Nodes();
               int nn_a = el.ShapeFunctions(dofPerNode_a);
               int nn_b = el.ShapeFunctions(dofPerNode_b);

               // get the element stiffness matrix
               int ndof_a = dofPerNode_a * nn_a;
               int ndof_b = dofPerNode_b * nn_b;
               ke.resize(ndof_a, ndof_b);

               // calculate element stiffness
               int nint = el.GaussPoints();

               // gauss weights
               double* w = el.GaussWeights();

               // nodal coordinates
               surf.GetNodalCoordinates(el, rt);

               // repeat over integration points
               ke.zero();
               for (int n = 0; n<nint; ++n)
               {
                  FESurfaceMaterialPoint& pt = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));

                  double* Gr = el.Gr(n);
                  double* Gs = el.Gs(n);

                  // tangents at integration point
                  pt.dxr = vec3d(0, 0, 0);
                  pt.dxs = vec3d(0, 0, 0);
                  for (int j = 0; j<neln; ++j)
                  {
                     pt.dxr += rt[j] * Gr[j];
                     pt.dxs += rt[j] * Gs[j];
                  }

                  // calculate stiffness component
                  for (int ia = 0; ia < nn_a; ++ia)
                  {
                     // shape function values
                     dof_a.index = ia;
                     dof_a.shape = el.H(order_a, n)[ia];
                     dof_a.shape_deriv_r = el.Gr(order_a, n)[ia];
                     dof_a.shape_deriv_s = el.Gs(order_a, n)[ia];

                     for (int jb = 0; jb < nn_b; ++jb)
                     {
                        // shape function values
                        dof_b.index = jb;
                        dof_b.shape = el.H(order_b, n)[jb];
                        dof_b.shape_deriv_r = el.Gr(order_b, n)[jb];
                        dof_b.shape_deriv_s = el.Gs(order_b, n)[jb];

                        // evaluate integrand
                        kab.zero();
                        sl->StiffnessIntegrand(pt, dof_a, dof_b, kab, factors);

                        // add it to the local element matrix
                        ke.adds(dofPerNode_a * ia, dofPerNode_b * jb, kab, -w[n]);
                     }
                  }
               }

               // assemble
               mfem::SparseMatrix& mfemK = *m_K->GetMfemMatrix();
               mfem::DenseMatrix mfem_elmat(*ke, vdofs_febio.Size(), vdofs_febio.Size());
               mfemK.AddSubMatrix(vdofs_febio, vdofs_febio, mfem_elmat, 0);
            }
         }
      }
   }

   // TODO: contact

   // Get MFEM version of FEBio stiffness matrix
   mfem::SparseMatrix *mfemK = m_K->GetMfemMatrix();
   if (!mfemK->Finalized())
   {
      mfemK->Finalize(0);
   }

   #ifdef MFEM_USE_MPI
   // Form parallel stiffness matrix
   mfem::OperatorHandle dK(mfem::Operator::Hypre_ParCSR);
   mfem::OperatorHandle Ph(mfem::Operator::Hypre_ParCSR);

   dK.MakeSquareBlockDiag(dfes->GetComm(), dfes->GlobalVSize(),dfes->GetDofOffsets(),mfemK);
   Ph.ConvertFrom(dfes->Dof_TrueDof_Matrix());

   pK.Clear();
   pK.MakePtAP(dK,Ph);
   mfem::HypreParMatrix *parmfemK = pK.As<mfem::HypreParMatrix>();
   #endif

   // if (disp_ess_tdof_list.Size() > 0)
   // {
   //    // eliminate essential displacement dofs
   //    #ifdef MFEM_USE_MPI
   //    parmfemK->EliminateBC(disp_ess_tdof_list,mfem::AbstractSparseMatrix::DIAG_ONE);
   //    #else
   //    MFEM_ABORT("TODO in serial");
   //    #endif
   // }

   sw_tan->Stop();

   //////////////////////////////////////////////////////////////////
   sw_fsitan->Start();

   #ifdef MFEM_USE_MPI
   // Get C matrix (fluid x solid dofs)
   #ifdef MFEM_USE_GSLIB
   mfem::HypreParMatrix *C = dynamic_cast<mfem::HypreParMatrix *>(&finder_fluid->GetCouplingMatrix());
   mfem::HypreParMatrix *Ct = dynamic_cast<mfem::HypreParMatrix *>(&finder_fluid->GetCouplingMatrixTranspose());
   #else
   mfem::HypreParMatrix *C = nullptr;
   mfem::HypreParMatrix *Ct = nullptr;
   #endif

   // C K Ct
   if (KCt) delete KCt;
   KCt = ParMult(parmfemK,Ct,true);
   if (CKCt) delete CKCt;
   CKCt = ParMult(C,KCt,true);

   sw_fsitan->Stop();
   return *CKCt;
   
   #else
   // Get C matrix (fluid x solid dofs)
   #ifdef MFEM_USE_GSLIB
   mfem::SparseMatrix *C = dynamic_cast<mfem::SparseMatrix *>(&finder_fluid->GetCouplingMatrix());
   mfem::SparseMatrix *Ct = dynamic_cast<mfem::SparseMatrix *>(&finder_fluid->GetCouplingMatrixTranspose());
   #else
   mfem::SparseMatrix *C = nullptr;
   mfem::SparseMatrix *Ct = nullptr;
   #endif

   // C K Ct
   if (KCt) delete KCt;
   KCt = mfem::Mult(*mfemK,*Ct);
   if (CKCt) delete CKCt;
   CKCt = mfem::Mult(*C,*KCt);

   sw_fsitan->Stop();
   return *CKCt;
   #endif

   // Eliminate essential dofs (TODO)
}

void ImmersedFEBioSolid::ReformStiffness()
{
   // // TODO: figure out why this doesnt work in parallel
   // if (reshape_stiffness)
   // {
   //    // setup sparse matrix profile the first time
   //    // update as needed (TODO: contact)
   //    m_K->CreateProfile(fem,dfes->GetVSize(),true);

   //    reshape_stiffness = false;
   // }

   // // // setup sparse matrix profile the first time
   // // // update as needed (TODO: contact)
   // // m_K->CreateProfile(fem,dfes->GetVSize(),true);
   
   // // zero the global matrix
   // m_K->Zero();



   if (reshape_stiffness)
   {
      // setup sparse matrix profile the first time
      // update as needed (TODO: contact)
      m_K->CreateProfile(fem,dfes->GetVSize(),true);

      reshape_stiffness = false;
   }

   mfem::SparseMatrix *mfemK = m_K->GetMfemMatrix();
   *mfemK = 0.0;
}

void ImmersedFEBioSolid::UpdateSolidStateFromRestart()
{
   ImmersedSolid::UpdateSolidStateFromRestart();
   UpdateFEBioDataFromRestart();
}

void ImmersedFEBioSolid::UpdateFEBioDataFromRestart()
{
   #ifndef MFEM_USE_MPI
   // // transfer submesh delvdot to global
   // delvdot_global_gf = 0.0;
   // submeshsolid_to_global_d->Transfer(temp_delvdot_ib_gf,delvdot_global_gf);
   // mfem::Vector &delvdot = delvdot_global_gf;

   // // transfer displacement to global
   // dnp1_gf.SetFromTrueDofs(dnp1);
   // *dnp1_global_gf = 0.0;
   // submeshsolid_to_global_d->Transfer(dnp1_gf,*dnp1_global_gf);
   // mfem::Vector &d = *dnp1_global_gf;

   MFEM_ABORT("TODO in serial");

   #else
   mfem::Vector &d = dn_gf;
   mfem::Vector &ddot = ddotn_gf;
   mfem::Vector &d2dot = d2dotn_gf;
   #endif
   
   // update nodal data
	FEMesh& mesh = fem->GetMesh();
   int dofVIMX = fem->GetDOFIndex("im_vx");
   int dofVIMY = fem->GetDOFIndex("im_vy");
   int dofVIMZ = fem->GetDOFIndex("im_vz");
   int NN = mesh.Nodes();
	for (int i=0; i<NN; ++i)
	{      
		FENode& ni = mesh.Node(i);
      
      int idx = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 0*NN : i*3 + 0;
      int idy = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 1*NN : i*3 + 1;
      int idz = dfes->GetOrdering() == mfem::Ordering::byNODES ? i + 2*NN : i*3 + 2;

      // acceleration
      ni.m_at(0) = d2dot[idx];
      ni.m_at(1) = d2dot[idy];
      ni.m_at(2) = d2dot[idz];

      // velocity
      vec3d vt = ni.get_vec3d(dofVIMX, dofVIMY, dofVIMZ);
      vt(0) = ddot[idx];
      vt(1) = ddot[idy];
      vt(2) = ddot[idz];
      ni.set_vec3d(dofVIMX, dofVIMY, dofVIMZ, vt);

      // position from displacement
      ni.m_rt(0) = d[idx] + ni.m_r0(0);
      ni.m_rt(1) = d[idy] + ni.m_r0(1);
      ni.m_rt(2) = d[idz] + ni.m_r0(2);
   }

   // update material point quantities
   // this should also be read in from the restart file
   // for (int i=0; i<mesh.Domains(); ++i)
   // {
   //    FEDomain& dom = mesh.Domain(i);
   //    MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
   //    if (mfd && mfd->IsActive())
   //    {
   //       FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(mfd->GetMaterial());
   //       assert(pme);
   //       for (int j=0; j<mfd->Elements(); ++j)
   //       {
   //          FESolidElement& el = mfd->Element(j);
   //          if (el.isActive())
   //          {
   //             int n = el.GaussPoints();
   //             for (int k = 0; k < n; ++k)
   //             {
   //                FEMaterialPoint& mp = *el.GetMaterialPoint(k);
   //                FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
   //                pt.m_Wt = pme->StrainEnergyDensity(mp);

   //                // TO DO: mp.Update(timeInfo);
   //             }
   //          }
   //       }
   //    }      
   // }
}

GridFunction *ImmersedFEBioSolid::ComputeNodallyAveragedVonMisesStress(GridFunction &d_gf)
{
   if (nodal_stress_fec == nullptr)
   {
      // create a linear fes for stress output
      int dim = pmesh->Dimension();
      nodal_stress_fec = new H1_FECollection(orderd, dim);
      #ifdef MFEM_USE_MPI
      nodal_stress_fes = new ParFiniteElementSpace(pmesh, nodal_stress_fec);
      #else
      nodal_stress_fes = new FiniteElementSpace(pmesh, nodal_stress_fec);
      #endif
      nodal_stress_gf.SetSpace(nodal_stress_fes);
   }
   nodal_stress_gf = 0.0;

   // AccumulateAndCountZones.
   Array<int> zones_per_vdof;
   zones_per_vdof.SetSize(nodal_stress_gf.FESpace()->GetVSize());
   zones_per_vdof = 0;

   FEMesh& mesh = fem->GetMesh();
   for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
      MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (mfd)
      {
         mfd->ComputeNodallyAveragedVonMisesStress(nodal_stress_gf, zones_per_vdof);
      }	
	}

   #ifdef MFEM_USE_MPI
   // average across processors
   // Count the zones globally
   auto par_nodal_stress_gf = dynamic_cast<ParGridFunction*>(&nodal_stress_gf);
   MFEM_ASSERT(par_nodal_stress_gf != nullptr, "Expected a ParGridFunction.");
   GroupCommunicator &gcomm = par_nodal_stress_gf->ParFESpace()->GroupComm();
   gcomm.Reduce<int>(zones_per_vdof, GroupCommunicator::Sum);
   gcomm.Bcast(zones_per_vdof);

   // Accumulate for all vdofs
   gcomm.Reduce<double>(nodal_stress_gf.GetData(), GroupCommunicator::Sum);
   gcomm.Bcast<double>(nodal_stress_gf.GetData());
   #endif

   // Compute means
   for (int i = 0; i < nodal_stress_gf.Size(); i++)
   {
      const int nz = zones_per_vdof[i];
      if (nz)
      {
         nodal_stress_gf(i) /= nz;
      }
   }

   return &nodal_stress_gf;
}

ImmersedFEBioSolid::~ImmersedFEBioSolid()
{
   delete m_R;
   delete d2dot_global_gf;
   delete dnp1_global_gf;
   delete m_K;

   for (int i=0; i<m_body_forces.size(); i++)
   {
      m_body_forces[i] = nullptr;
   }

   #ifdef MFEM_USE_MPI
   pK.Clear();
   #endif
}