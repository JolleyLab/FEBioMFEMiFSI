#include "MFEMImmersedFEBioSolid.h"
#include <FECore/FEDomain.h>
#include <FECore/FEModelLoad.h>
#include <FEBioMech/FEElasticDomain.h>
#include "MFEMImmersedElasticDomain.h"
#include "mfem.hpp"
#include <general/forall.hpp>

using namespace mfem;
using namespace immersed_solid_serial;

ImmersedFEBioSolid::ImmersedFEBioSolid(FEModel* fem_, mfem::Mesh *solid_mesh, int order_, mfem::Coefficient *init_rho_s_, double fluid_density_, double mu_f_, std::vector<AccelTerm_T> *accel_,
   mfem::FiniteElementSpace *vfes_, mfem::Mesh *fluid_mesh, MaterialType *mat_, int orderingtype_, int form_type_, bool subtract_visc_)
: ImmersedSolid(solid_mesh, order_, init_rho_s_, fluid_density_, mu_f_, accel_, vfes_, fluid_mesh, mat_, orderingtype_, form_type_, subtract_visc_), fem(fem_)
{
   int dfes_truevsize = dfes->GetTrueVSize();

   // setup and initialize additional variables
   dnp1_gf.SetSpace(dfes);
   dnp1_gf = 0.0;
   ddotnp1_gf.SetSpace(dfes);
   ddotnp1_gf = 0.0;
   d2dotnp1_gf.SetSpace(dfes);
   d2dotnp1_gf = 0.0;
   
   // setup residual vector
   R.resize(dfes_truevsize,0.0);
   m_R = new MFEMGlobalVector(fem,dfes_truevsize,R,Fr);

   // setup linear system to assemble FEBio into
   CRSSparseMatrix *pS = new CRSSparseMatrix(0);
   if (m_K) delete m_K;
   m_K = new MFEMGlobalMatrix(pS); // ps owned by MFEMGlobalMatrix
}


void ImmersedFEBioSolid::Setup(double dt, double alpha_m_, double alpha_f_, double gamma_, double beta_)
{
   mfem::out << "Setup Immersed Solid" << std::endl;
   PrintInfo();

   this->dt = dt;
   
   // setup point search
   #ifdef MFEM_USE_GSLIB
   finder_fluid = new ImmersedFSI_FindPts();
   finder_fluid->Setup(*fluid_pmesh);
   #endif

   true_interp_points_in_solid.SetSize(dfes->GetTrueVSize());
   true_interp_points_in_solid = 0.0;

   interp_points_in_solid = *pmesh->GetNodes();
   mfem::GridFunction *nodes = pmesh->GetNodes();
   nodes->GetTrueDofs(true_interp_points_in_solid);

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
   FEMesh &mesh = fem->GetMesh();
   for (int nd = 0; nd<mesh.Domains(); ++nd)
   {
      FEDomain& dom = mesh.Domain(nd);
      MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (mfd)
      {
         mfd->SetFESpace(dfes);
         mfd->SetFluidProperties(fluid_density,mu_f);
         mfd->SetAdditionalData(subtract_fluid_visc,dynamic_cast<mfem::SubMesh*>(pmesh),ddot_global_gf->FESpace());
      }			
   }

   // setup VMS stab scaling
   if (stabs_scaling)
   {
      // setup interpolation points (only boundary dofs of the solid)
      mfem::Array<int> bdr_attr(pmesh->bdr_attributes.Max());
      bdr_attr = 1;
      dfes->GetEssentialTrueDofs(bdr_attr,solid_bdr_tdofs);
      solid_bdr_nodes.SetSize(solid_bdr_tdofs.Size());

      const int Ndofs = solid_bdr_tdofs.Size();
      const auto id = solid_bdr_tdofs.Read();
      auto FV_nodes = true_interp_points_in_solid.Read();
      auto FV = solid_bdr_nodes.ReadWrite();
      mfem::forall(Ndofs, [=] MFEM_HOST_DEVICE (int i) { FV[i] = FV_nodes[id[i]]; });
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
   ddotnp1_gf.GetTrueDofs(ddotnp1);

   // d2dotnp1 = vdot_np1
   temp_del_gf.SetFromTrueDofs(vdotf);
   #ifdef MFEM_USE_GSLIB
   finder_fluid->Interpolate(interp_points_in_solid,temp_del_gf,d2dotnp1_gf,orderingtype);
   #endif
   // d2dotnp1_gf.GetTrueDofs(d2dotnp1);

   // transfer gf defined on solid submesh to global mesh (FEBio representation)
   *ddot_global_gf = 0.0; // not sure if needed?
   *d2dot_global_gf = 0.0;
   submeshsolid_to_global_d->Transfer(ddotnp1_gf,*ddot_global_gf);
   submeshsolid_to_global_d->Transfer(d2dotnp1_gf,*d2dot_global_gf);

   int dofVIMX = fem->GetDOFIndex("im_vx");
   int dofVIMY = fem->GetDOFIndex("im_vy");
   int dofVIMZ = fem->GetDOFIndex("im_vz");
   mfem::Vector &ddot = *ddot_global_gf;
   mfem::Vector &d2dot = *d2dot_global_gf;
   mfem::Vector &d = *dnp1_global_gf;

   // Update nodal data, which stores values at n and np1 (from SolidSolver2::PrepStep())   
   // switch this to a loop over the elements in the immersed domain and only update the nodal
   // values in those elements
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

   // Transfer global dnp1 to submesh
   dnp1_gf = 0.0; // not sure if needed?
   global_to_submeshsolid_d->Transfer(*dnp1_global_gf,dnp1_gf);
   dnp1_gf.GetTrueDofs(dnp1);

   // // Open a file to write the vector
   // std::ofstream outfile("dnp1_global.txt", std::ios::out);

   // int precision_out = 6;
   // int width_out = precision_out + 6;

   // // Print the matrix ke to the file
   // for (int i = 0; i < dnp1.Size(); ++i) 
   // {
   //    outfile << std::setw(width_out) << std::setprecision(precision_out) << dnp1(i);
   //    outfile << std::endl;
   // }
   // outfile << std::endl; // Add an extra line between matrices
   // // Close the file
   // outfile.close();

   // // int dofVIMX = fem->GetDOFIndex("im_vx");
   // // int dofVIMY = fem->GetDOFIndex("im_vy");
   // // int dofVIMZ = fem->GetDOFIndex("im_vz");
   // // Vector &ddot = ddotnp1_gf;
   // // Vector &d2dot = d2dotnp1_gf;
   // // Vector &d = dnp1_gf;

   // dnp1_gf = 0.0;

   // // FETimeInfo& tp = fem->GetTime();
	// // FEMesh& mesh = fem->GetMesh();
   // Array<int> vdofs;
   // for (int i=0; i<mesh.Domains(); ++i) 
	// {
	// 	FEDomain& dom = mesh.Domain(i);
   //    MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
   //    if (mfd)
   //    {
   //       SubMesh &solid_submesh = *mfd->GetSubMesh();
   //       for (int e=0; e<mfd->Elements(); ++e)
   //       {
   //          FESolidElement& el = mfd->Element(e);
   //          if (el.isActive())
   //          {
   //             int submesh_id = solid_submesh.GetSubMeshElementFromParent(el.GetID()-1);
   //             dfes->GetElementVDofs(submesh_id,vdofs);
   //             int nd = el.Nodes();

   //             for (int n=0; n<el.Nodes(); ++n)
   //             {
   //                FENode& node = mesh.Node(el.m_node[n]);
                  
   //                // store prev state
   //                node.m_rp = node.m_rt;
   //                node.m_vp = node.get_vec3d(dofVIMX, dofVIMY, dofVIMZ);
   //                node.m_ap = node.m_at;
   //                node.UpdateValues();

   //                // solid predictors
   //                int idx = vdofs[n];
   //                int idy = vdofs[n+nd];
   //                int idz = vdofs[n+2*nd];

   //                // acceleration
   //                node.m_at(0) = d2dotnp1_gf[idx];
   //                node.m_at(1) = d2dotnp1_gf[idy];
   //                node.m_at(2) = d2dotnp1_gf[idz];

   //                // velocity
   //                vec3d vs(ddotnp1_gf[idx],ddotnp1_gf[idy],ddotnp1_gf[idz]);
   //                node.set_vec3d(dofVIMX, dofVIMY, dofVIMZ, vs);

   //                // position from displacement predictor
   //                double fac_n = 0.5*dt*dt*(1.0 - 2.0*beta);
   //                double fac_np1 = 0.5*dt*dt*(2.0*beta);
   //                vec3d dn(node.m_rt(0), node.m_rt(1), node.m_rt(2));
   //                dn -= node.m_r0; // compute dn = xt - x0
   //                double dnp1_x = dn(0) + dt*node.m_vp(0) + fac_n*node.m_ap(0) + fac_np1*node.m_at(0);
   //                double dnp1_y = dn(1) + dt*node.m_vp(1) + fac_n*node.m_ap(1) + fac_np1*node.m_at(1);
   //                double dnp1_z = dn(2) + dt*node.m_vp(2) + fac_n*node.m_ap(2) + fac_np1*node.m_at(2);
   //                node.m_rt(0) = dnp1_x + node.m_r0(0);
   //                node.m_rt(1) = dnp1_y + node.m_r0(1);
   //                node.m_rt(2) = dnp1_z + node.m_r0(2);

   //                // set mfem dnp1
   //                dnp1_gf[idx] = dnp1_x;
   //                dnp1_gf[idy] = dnp1_y;
   //                dnp1_gf[idz] = dnp1_z;
   //             }
   //          }
   //       }
   //    }	
	// }
   // dnp1_gf.GetTrueDofs(dnp1);
   
   // // Open a file to write the vector
   // std::ofstream outfile2("dnp1_local.txt", std::ios::out);

   // // int precision_out = 6;
   // // int width_out = precision_out + 6;

   // // Print the matrix ke to the file
   // for (int i = 0; i < dnp1.Size(); ++i) 
   // {
   //    outfile2 << std::setw(width_out) << std::setprecision(precision_out) << dnp1(i);
   //    outfile2 << std::endl;
   // }
   // outfile2 << std::endl; // Add an extra line between matrices
   // // Close the file
   // outfile2.close();

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
}

void ImmersedFEBioSolid::Initiator()
{
   // update element data (material point data, which is data stored
   // at the integration points) to the generalized midpoint
	fem->Update();

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
   UpdatePointSearch(dtemp_gf,dg);

   if (stabs_scaling)
   {
      #ifdef MFEM_USE_GSLIB
      MFEM_ASSERT(dfes->GetOrdering() == s_gf_fluid->FESpace()->GetOrdering(),"Ordering of solid and fluid gridfunctions must match for special VMS stabs scaling");
      
      finder_fluid->TransferIndicator(solid_bdr_nodes,*s_gf_fluid,dfes->GetOrdering());
      // note internally s_gf_fluid is zeroed prior to interpolation of s_gf to its fluid counterpart
      #endif
   }
}

void ImmersedFEBioSolid::Correctors(mfem::Vector &del)
{
   // project delvdot to solid mesh and update data
   temp_del_gf.SetFromTrueDofs(del);
   #ifdef MFEM_USE_GSLIB
   finder_fluid->Interpolate(interp_points_in_solid,temp_del_gf,temp_delvdot_ib_gf,orderingtype);
   #endif
   temp_delvdot_ib_gf.GetTrueDofs(delvdot_ib);

   // d2dotnp1.Add(-1.0,delvdot_ib);
   ddotnp1.Add(-gamma*dt,delvdot_ib); // only updated to keep ddotn_gf current
   dnp1.Add(-beta*dt*dt,delvdot_ib); // for moving mfem mesh

   ////////////////////////////////////////////////////////////////

   // transfer submesh delvdot to global.
   delvdot_global_gf = 0.0;
   submeshsolid_to_global_d->Transfer(temp_delvdot_ib_gf,delvdot_global_gf);
   mfem::Vector &delvdot = delvdot_global_gf;

   // transfer displacement to global
   dnp1_gf.SetFromTrueDofs(dnp1);
   *dnp1_global_gf = 0.0;
   submeshsolid_to_global_d->Transfer(dnp1_gf,*dnp1_global_gf);
   mfem::Vector &d = *dnp1_global_gf;
   
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

      // if (i == 3626 || i == 3643 || i == 3783 || i == 3662 || i == 4083 || i == 4100 || i == 4240 || i == 4119)
      // {
      //    mfem::out << "--------------------------------\n";
      //    mfem::out << "Global node id = " << i << std::endl;
      //    mfem::out << "Node id = " << i << std::endl;
      //    mfem::out << "ni.m_at[0] = " << ni.m_at(0) << std::endl;
      //    mfem::out << "ni.m_at[1] = " << ni.m_at(1) << std::endl;
      //    mfem::out << "ni.m_at[2] = " << ni.m_at(2) << std::endl;
      //    mfem::out << "ni.m_rt[0] = " << ni.m_rt(0) << std::endl;
      //    mfem::out << "ni.m_rt[1] = " << ni.m_rt(1) << std::endl;
      //    mfem::out << "ni.m_rt[2] = " << ni.m_rt(2) << std::endl;
      //    mfem::out << "ni.m_vt[0] = " << ni.get(dofVIMX) << std::endl;
      //    mfem::out << "ni.m_vt[1] = " << ni.get(dofVIMY) << std::endl;
      //    mfem::out << "ni.m_vt[2] = " << ni.get(dofVIMZ) << std::endl;
      // }
   }

   // int dofVIMX = fem->GetDOFIndex("im_vx");
   // int dofVIMY = fem->GetDOFIndex("im_vy");
   // int dofVIMZ = fem->GetDOFIndex("im_vz");
   // Vector &delvdot_loc = temp_delvdot_ib_gf;
   // dnp1_gf.SetFromTrueDofs(dnp1);
   // Vector &d_loc = dnp1_gf;

   // FETimeInfo& tp = fem->GetTime();
	// FEMesh& mesh = fem->GetMesh();
   // Array<int> vdofs;
   // for (int i=0; i<mesh.Domains(); ++i) 
	// {
	// 	FEDomain& dom = mesh.Domain(i);
   //    MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
   //    if (mfd)
   //    {
   //       SubMesh &solid_submesh = *mfd->GetSubMesh();
   //       for (int e=0; e<mfd->Elements(); ++e)
   //       {
   //          FESolidElement& el = mfd->Element(e);
   //          if (el.isActive())
   //          {
   //             int submesh_id = solid_submesh.GetSubMeshElementFromParent(el.GetID()-1);
   //             dfes->GetElementVDofs(submesh_id,vdofs);
   //             int nd = el.Nodes();

   //             // mfem::out << "Element id = " << submesh_id << std::endl;
   //             // auto curr_elem = dfes->GetFE(submesh_id);
   //             // auto trans = dfes->GetElementTransformation(submesh_id);
   //             // auto Nodes = curr_elem->GetNodes();
   //             // mfem::Vector curr_pos(3);
   //             // for(int j=0; j<Nodes.GetNPoints(); ++j)
   //             // {
   //             //    auto ip_node = Nodes.IntPoint(j);
   //             //    trans->SetIntPoint(&ip_node);
   //             //    trans->Transform(ip_node, curr_pos);
   //             //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
   //             // }
   //             // mfem::out << "------------------------\n";

   //             for (int n=0; n<el.Nodes(); ++n)
   //             {
   //                int glob_nodeid = el.m_node[n];
   //                FENode& node = mesh.Node(el.m_node[n]);

   //                // mfem::out << std::setprecision(8) << "CurrNode[" << n << "]=(" << node.m_rt(0) << "," << node.m_rt(1) << "," << node.m_rt(2) << ")" << std::endl;

   //                int idx = vdofs[n];
   //                int idy = vdofs[n+nd];
   //                int idz = vdofs[n+2*nd];

   //                // acceleration
   //                node.m_at(0) -= delvdot_loc[idx];
   //                node.m_at(1) -= delvdot_loc[idy];
   //                node.m_at(2) -= delvdot_loc[idz];

   //                // velocity
   //                vec3d vt = node.get_vec3d(dofVIMX, dofVIMY, dofVIMZ);
   //                vt(0) -= gamma*dt*delvdot_loc[idx];
   //                vt(1) -= gamma*dt*delvdot_loc[idy];
   //                vt(2) -= gamma*dt*delvdot_loc[idz];
   //                node.set_vec3d(dofVIMX, dofVIMY, dofVIMZ, vt);

   //                // position from updated displacement
   //                node.m_rt(0) = d_loc[idx] + node.m_r0(0);
   //                node.m_rt(1) = d_loc[idy] + node.m_r0(1);
   //                node.m_rt(2) = d_loc[idz] + node.m_r0(2);

   //                // print node data
                  // mfem::out << "--------------------------------\n";
                  // mfem::out << "Global node id = " << glob_nodeid << std::endl;
   //                mfem::out << "node.m_at[0] = " << node.m_at(0) << std::endl;
   //                mfem::out << "node.m_at[1] = " << node.m_at(1) << std::endl;
   //                mfem::out << "node.m_at[2] = " << node.m_at(2) << std::endl;
   //                mfem::out << "node.m_rt[0] = " << node.m_rt(0) << std::endl;
   //                mfem::out << "node.m_rt[1] = " << node.m_rt(1) << std::endl;
   //                mfem::out << "node.m_rt[2] = " << node.m_rt(2) << std::endl;
   //                mfem::out << "node.m_vt[0] = " << node.get(dofVIMX) << std::endl;
   //                mfem::out << "node.m_vt[1] = " << node.get(dofVIMY) << std::endl;
   //                mfem::out << "node.m_vt[2] = " << node.get(dofVIMZ) << std::endl;
   //             }
   //          }
   //       }
   //    }	
	// }
}

void ImmersedFEBioSolid::Update()
{
   // update dn_gf and ddotn_gf in mfem,
   dn = dnp1;
   dn_gf.SetFromTrueDofs(dn);

   // update ddotn_gf in mfem
   ddotn = ddotnp1;
   ddotn_gf.SetFromTrueDofs(ddotn);

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

	// calculate residual (note: no surface forces considered)
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

   // TODO: contact

   sw_res->Stop();

   ////////////////////////////////////////////////////////
   // Get pointer to mfem::Vector version of m_R
   mfem::Vector &residual_vector = m_R->GetMfemVector();

   sw_fsires->Start();
   // Spread to fluid nodes
   #ifdef MFEM_USE_GSLIB
   mfem::SparseMatrix *C = dynamic_cast<mfem::SparseMatrix *>(&finder_fluid->GetCouplingMatrix());
   #else
   mfem::SparseMatrix *C = nullptr;
   #endif
   C->Mult(residual_vector,lambda_fsi_vec);   
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
   for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
      MFEMImmersedElasticDomain *mfd = dynamic_cast<MFEMImmersedElasticDomain*>(&dom);
      if (mfd)
      {
         mfd->AssembleDomainStiffness(*m_K,am,av,af,m_body_forces);
      }	
	}

   // TODO: contact

   sw_tan->Stop();
   //////////////////////////////////////////////////////////////////
   sw_fsitan->Start();

   // Get MFEM version of FEBio stiffness matrix
   mfem::SparseMatrix *mfemK = m_K->GetMfemMatrix();

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

   // Eliminate essential dofs (TODO)
}

void ImmersedFEBioSolid::ReformStiffness()
{
   if (reshape_stiffness)
   {
      // setup sparse matrix profile the first time
      // update as needed (TODO: contact)
      m_K->CreateProfile(fem,dfes->GetTrueVSize(),true);

      reshape_stiffness = false;
   }
   
   // zero the global matrix
   m_K->Zero();
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
}