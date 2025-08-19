#include "MFEMImmersedElasticDomain.h"
#include <FEBioMech/FEBodyForce.h>
#include <FECore/FEModel.h>
#include "MFEMGlobalTypes.h"

#include <fstream>
#include <iomanip>

MFEMImmersedElasticDomain::MFEMImmersedElasticDomain(FEModel* pfem) : FESolidDomain(pfem), m_dof(pfem)
{
   fluid_density = 0.0;
   fluid_dyn_viscosity = 0.0;
}

MFEMImmersedElasticDomain::~MFEMImmersedElasticDomain() noexcept
{
   
}

void MFEMImmersedElasticDomain::SetFluidProperties(double rho, double mu)
{
   fluid_density = rho;
   fluid_dyn_viscosity = mu;
}

#ifdef MFEM_USE_MPI
void MFEMImmersedElasticDomain::SetAdditionalData(bool subtract_visc_fluid_)
{
   auto tp = GetFEModel()->GetTime();
   dt = tp.timeIncrement;
   m_gamma = tp.gamma;

   subtract_visc_fluid = subtract_visc_fluid_;

   solid_submesh = nullptr;
   gfes = nullptr;
}
#else
void MFEMImmersedElasticDomain::SetAdditionalData(bool subtract_visc_fluid_, mfem::SubMesh* solid_submesh_, mfem::FiniteElementSpace *gfes_)
{
   auto tp = GetFEModel()->GetTime();
   dt = tp.timeIncrement;
   m_gamma = tp.gamma;

   subtract_visc_fluid = subtract_visc_fluid_;

   solid_submesh = solid_submesh_;
   gfes = gfes_;
}
#endif

void MFEMImmersedElasticDomain::Activate()
{  
   FEModel& fem = *GetFEModel();
   int dofUX = fem.GetDOFIndex("im_ux");
   int dofUY = fem.GetDOFIndex("im_uy");
   int dofUZ = fem.GetDOFIndex("im_uz");
	for (int i=0; i<Nodes(); ++i)
	{
		FENode& node = Node(i);
		if (node.HasFlags(FENode::EXCLUDE) == false)
		{
			if (node.m_rid < 0)
			{
				node.set_active(dofUX);
				node.set_active(dofUY);
				node.set_active(dofUZ);
			}
		}
	}
}

//! Initialize element data
void MFEMImmersedElasticDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    m_alphaf = timeInfo.alphaf;
    m_alpham = timeInfo.alpham;
    m_beta = timeInfo.beta;

// #pragma omp parallel for
	for (int i=0; i<Elements(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		if (el.isActive())
		{
			int n = el.GaussPoints();
			for (int j = 0; j < n; ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				pt.m_Wp = pt.m_Wt;

				mp.Update(timeInfo);
			}
		}
	}
}

void MFEMImmersedElasticDomain::Update(const FETimeInfo& tp)
{
	bool berr = false;
	int NE = Elements();
	// #pragma omp parallel for shared(NE, berr)
	for (int i=0; i<NE; ++i)
	{
		try
		{
			FESolidElement& el = Element(i);
			if (el.isActive())
			{
				UpdateElementStress(i, tp);
			}
		}
		catch (NegativeJacobian e)
		{
			// #pragma omp critical
			{
				// reset the logfile mode
				berr = true;
				if (e.DoOutput()) feLogError(e.what());
			}
		}
	}

	if (berr) throw NegativeJacobianDetected();
}

//! Update element state data (mostly stresses, but some other stuff as well)
//! \todo Remove the remodeling solid stuff
void MFEMImmersedElasticDomain::UpdateElementStress(int iel, const FETimeInfo& tp)
{
   double dt =tp.timeIncrement;
    
	// get the solid element
	FESolidElement& el = m_Elem[iel];

	// get the number of integration points
	int nint = el.GaussPoints();

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
   const int NELN = FEElement::MAX_NODES;
   vec3d r[NELN], v[NELN], a[NELN];
	GetCurrentNodalCoordinates(el, r, m_alphaf);

	// update dynamic quantities
	if (m_update_dynamic)
	{
      int dofVIMX = GetFEModel()->GetDOFIndex("im_vx");
      int dofVIMY = GetFEModel()->GetDOFIndex("im_vy");
      int dofVIMZ = GetFEModel()->GetDOFIndex("im_vz");

		for (int j = 0; j<neln; ++j)
		{
			FENode& node = m_pMesh->Node(el.m_node[j]);
			v[j] = node.get_vec3d(dofVIMX, dofVIMY, dofVIMZ)*m_alphaf + node.m_vp*(1 - m_alphaf);
			a[j] = node.m_at*m_alpham + node.m_ap*(1 - m_alpham);
		}
	}

	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// material point coordinates
		mp.m_rt = el.Evaluate(r, n); // sets current position at int point to midpoint value

		// get the deformation gradient and determinant at intermediate time
      double Jt;
      mat3d Ft, Fp;
      Jt = defgrad(el, Ft, n);
      defgradp(el, Fp, n);
		if (m_alphaf == 1.0)
		{
			pt.m_F = Ft;
         pt.m_J = Jt;
		}
		else
		{
			pt.m_F = Ft*m_alphaf + Fp*(1-m_alphaf);
         pt.m_J = pt.m_F.det();
		}

      // mat3d Fi = pt.m_F.inverse();
      // pt.m_L = (Ft - Fp)*Fi / dt;
      pt.m_L = gradient(el, v, n, m_alphaf); // gradv at int point at the midpoint
		if (m_update_dynamic)
		{
         // sets current velocity and acceleration at int point to midpoints
			pt.m_v = el.Evaluate(v, n);
			pt.m_a = el.Evaluate(a, n);
		}

      // update specialized material points
      m_pMat->UpdateSpecializedMaterialPoints(mp, tp);
        
		// calculate the stress at this material point
		pt.m_s = (m_secant_stress ? m_pMat->SecantStress(mp) : m_pMat->Stress(mp));
        
      // adjust stress for strain energy conservation
		// (Apply only for mid-point rule)
		if (m_alphaf == 0.5)
		{
         MFEM_ABORT("TO DO");
			FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(m_pMat);

			// evaluate strain energy at current time
			mat3d Ftmp = pt.m_F;
			double Jtmp = pt.m_J;
			pt.m_F = Ft;
			pt.m_J = Jt;
			pt.m_Wt = pme->StrainEnergyDensity(mp);
			pt.m_F = Ftmp;
			pt.m_J = Jtmp;

         mat3ds D = pt.RateOfDeformation();
         double D2 = D.dotdot(D);
			if (D2 > std::numeric_limits<double>::epsilon())
			{
				pt.m_s += D * (((pt.m_Wt - pt.m_Wp) / (dt * pt.m_J) - pt.m_s.dotdot(D)) / D2);
			}
      }
   }
}

void MFEMImmersedElasticDomain::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*3);

   int dofUX = GetFEModel()->GetDOFIndex("im_ux");
   int dofUY = GetFEModel()->GetDOFIndex("im_uy");
   int dofUZ = GetFEModel()->GetDOFIndex("im_uz");
	for (int i=0; i<N; ++i)
	{
		FENode& node = m_pMesh->Node(el.m_node[i]);
		vector<int>& id = node.m_ID;
      
      lm[3*i  ] = id[dofUX];
		lm[3*i+1] = id[dofUY];
		lm[3*i+2] = id[dofUZ];
	}
}

mfem::Array<int> &MFEMImmersedElasticDomain::TransformVDofs(mfem::Array<int> &vdofs, int nd)
{
   febio_vdofs.SetSize(vdofs.Size());
   for (int i=0; i < nd; ++i)
   {
      febio_vdofs[3*i]   = vdofs[i];
      febio_vdofs[3*i+1] = vdofs[i+nd];
      febio_vdofs[3*i+2] = vdofs[i+2*nd];
   }
   return febio_vdofs;
}

//-----------------------------------------------------------------------------
// RESIDUAL
void MFEMImmersedElasticDomain::AssembleDomainResidual(MFEMGlobalVector& R, std::vector<MFEMBodyForce*> &bfs)
{
   // Loop over elements
	int NE = Elements();
	for (int i=0; i<NE; ++i)
	{
		// get the element
		FESolidElement& el = m_Elem[i];

      // get submesh element from parent
      // febio element id is 1-based, while mfem element id is 0-based
      // so we need to subtract 1 when getting the mfem element id
      
      int submesh_id = -1;
      if (solid_submesh)
      {
         // serial
         submesh_id = solid_submesh->GetSubMeshElementFromParent(el.GetID()-1);
      }         
      else
      {
         // parallel
         submesh_id = el.GetID()-1;
         // submesh_id = el.GetID();
      }

      // const mfem::Array<int>& parent_elem_ids = solid_submesh->GetParentElementIDMap();
      // int global_id = parent_elem_ids[submesh_id];

      // // auto parent_mesh = solid_submesh->GetParent();
      // // const_cast<mfem::Mesh*>(parent_mesh)->EnsureNodes();
      // // auto mesh_space = parent_mesh->GetNodes()->FESpace();

      // // print out nodal coordinates of each node on the mfem mesh
      // auto Nodes = gfes->GetFE(global_id)->GetNodes();
      // auto eltrans = gfes->GetElementTransformation(global_id);
      // mfem::out << "MFEM Nodes = " << "\n";
      // mfem::Vector curr_pos(3);
      // for(int j=0; j<Nodes.GetNPoints(); ++j)
      // {
      //    auto ip_node = Nodes.IntPoint(j);
      //    eltrans->SetIntPoint(&ip_node);
      //    eltrans->Transform(ip_node, curr_pos);
      //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
      // }

      // mfem::Array<int> elem_vert;
      // gfes->GetElementVertices(global_id, elem_vert);
      // elem_vert.Print();
      // mfem::out << "MFEM Nodes = " << "\n";
      // // mfem::Vector curr_pos(3);
      // for(int j=0; j<elem_vert.Size(); ++j)
      // {
      //    double* node = gfes->GetMesh()->GetVertex(elem_vert[j]);
      //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << node[0] << "," << node[1] << "," << node[2] << ")" << std::endl;
      // }

      // // print out nodal coordinates of each node on the mfem mesh
      // auto Nodes_mesh = mesh_space->GetFE(global_id)->GetNodes();
      // auto eltrans_mesh = mesh_space->GetElementTransformation(global_id);
      // mfem::out << "MFEM Mesh Nodes = " << "\n";
      // for(int j=0; j<Nodes_mesh.GetNPoints(); ++j)
      // {
      //    auto ip_node = Nodes_mesh.IntPoint(j);
      //    eltrans_mesh->SetIntPoint(&ip_node);
      //    eltrans_mesh->Transform(ip_node, curr_pos);
      //    mfem::out << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
      // }

      // mfem::out << "Element id = " << submesh_id << std::endl;
      // auto curr_elem = dfes->GetFE(submesh_id);
      // auto trans = dfes->GetElementTransformation(submesh_id);
      // auto Nodes = curr_elem->GetNodes();
      // mfem::Vector curr_pos(3);
      // for(int j=0; j<Nodes.GetNPoints(); ++j)
      // {
      //    auto ip_node = Nodes.IntPoint(j);
      //    trans->SetIntPoint(&ip_node);
      //    trans->Transform(ip_node, curr_pos);
      //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
      // }

		if (el.isActive()) {
			// element force vector
			vector<double> fe;
         // std::vector<int> lm;
			mfem::Array<int> vdofs;

			// get the element force vector and initialize it to zero
			int ndof = 3 * el.Nodes();
			fe.assign(ndof, 0);

         // internal force 
			ElementInternalForce(el, fe);

         // inertial force (acceleration - body force)
         ElementInertialForce(el, fe, bfs);

         if (subtract_visc_fluid)
         {
            // artificial fluid viscous force
            ElementArtificialFluidForce(el, fe);
         }

         // // get the element's LM vector
			// UnpackLM(el, lm);
         // mfem::out << "\nLM = ";
         // for (int ii=0; ii<lm.size(); ++ii)
         // {
         //    mfem::out << lm[ii] << " ";
         // }
         // mfem::out << std::endl;

			// assemble element 'fe'-vector into global R vector
         dfes->GetElementVDofs(submesh_id,vdofs);
         mfem::Array<int> vdofs_febio = TransformVDofs(vdofs,el.Nodes()); // reorder to match febio elem vector [xyz,xyz,xyz,...]

         // // mfem::out << "\nVDOFS = ";
         // // for (int ii=0; ii<vdofs_febio.Size(); ++ii)
         // // {
         // //    mfem::out << vdofs_febio[ii] << " ";
         // // }
         // // mfem::out << std::endl;

         // mfem::out << "Element Vector = " << std::endl;
         // mfem::Vector fe_temp(fe.data(),fe.size());
         // for (int ii=0; ii<el.Nodes(); ++ii)
         // {
         //    mfem::out << fe_temp(3*ii) << " ";
         //    mfem::out << fe_temp(3*ii + 1) << " ";
         //    mfem::out << fe_temp(3*ii + 2) << " ";
         //    mfem::out << std::endl;
         // }
         // // mfem::Vector& fe_mfem = TransformElemVec(fe,ndof);
         // // mfem::out << "Element Vector MFEM = " << std::endl;
         // // for (int ii=0; ii<el.Nodes(); ++ii)
         // // {
         // //    mfem::out << fe_mfem(ii) << " ";
         // //    mfem::out << fe_mfem(ii + ndof) << " ";
         // //    mfem::out << fe_mfem(ii + 2*ndof) << " ";
         // //    mfem::out << std::endl;
         // // }

         // // mfem::Vector& fe_mfem = TransformElemVec(fe,el.Nodes());

         // // mfem::out << "VDOFS = " << std::endl;
         // // vdofs.Print();
         // // for (int ii=0; ii<el.Nodes(); ++ii)
         // // {
         // //    mfem::out << "Node " << ii << " = (" << vdofs[ii] << ", " << vdofs[ii + el.Nodes()] << ", " << vdofs[ii + 2*el.Nodes()] << ")" << std::endl;
         // // }

         // // mfem::out << "VDOFS_FEBIO = " << std::endl;
         // // vdofs_febio.Print();
         // // for (int ii=0; ii<el.Nodes(); ++ii)
         // // {
         // //    mfem::out << "Node " << ii << " = (" << vdofs_febio[3*ii] << ", " << vdofs_febio[3*ii + 1] << ", " << vdofs_febio[3*ii + 2] << ")" << std::endl;
         // // }

			R.Assemble(vdofs_febio,fe);
		}
	}
}

void MFEMImmersedElasticDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3];

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		double detJt = (m_update_dynamic ? invjact(el, Ji, n, m_alphaf) : invjact(el, Ji, n));

		detJt *= gw[n];

		// get the stress vector for this integration point
        const mat3ds& s = pt.m_s;

		const double* Gr = el.Gr(n);
		const double* Gs = el.Gs(n);
		const double* Gt = el.Gt(n);

		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			double Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			double Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			double Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] += ( Gx*s.xx() +
				           Gy*s.xy() +
					       Gz*s.xz() )*detJt;

			fe[3*i+1] += ( Gy*s.yy() +
				           Gx*s.xy() +
					       Gz*s.yz() )*detJt;

			fe[3*i+2] += ( Gz*s.zz() +
				           Gy*s.yz() +
					       Gx*s.xz() )*detJt;
		}
	}
}

void MFEMImmersedElasticDomain::ElementInertialForce(FESolidElement& el, vector<double>& fe, std::vector<MFEMBodyForce*>& bfs)
{
   int nint = el.GaussPoints();
   int neln = el.Nodes();
   
   double*    gw = el.GaussWeights();

   vec3d inertial_force;

   // repeat for all integration points
   for (int n=0; n<nint; ++n)
   {
      FEMaterialPoint& mp = *el.GetMaterialPoint(n);
      FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
      double dens = m_pMat->Density(mp) - fluid_density*pt.m_J;
      double J0 = detJ0(el, n)*gw[n];

      // inertial force
      // acceleration
      inertial_force(0) = pt.m_a(0);
      inertial_force(1) = pt.m_a(1);
      inertial_force(2) = pt.m_a(2);
      // body forces
      for (auto &bf : bfs)
      {
         // get the force
         inertial_force(0) -= bf->m_bx(mp);
         inertial_force(1) -= bf->m_by(mp);
         inertial_force(2) -= bf->m_bz(mp);
      }
      
      double* H = el.H(n);
      for (int i=0; i<neln; ++i)
      {
         double tmp = H[i]*J0*dens;
         fe[3*i  ] += tmp*inertial_force.x;
         fe[3*i+1] += tmp*inertial_force.y;
         fe[3*i+2] += tmp*inertial_force.z;
      }
   }
}

void MFEMImmersedElasticDomain::ElementArtificialFluidForce(FESolidElement& el, vector<double>& fe)
{
	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3];

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

   // // velocity gradient
   // // TODO, create immersed material point class that computes gradv and viscous stress directly!
   // vec3d v[FEElement::MAX_NODES];
   // int dofVIMX = GetFEModel()->GetDOFIndex("im_vx");
   // int dofVIMY = GetFEModel()->GetDOFIndex("im_vy");
   // int dofVIMZ = GetFEModel()->GetDOFIndex("im_vz");
   // for (int j=0; j<neln; ++j) 
   // {
   //    FENode& node = m_pMesh->Node(el.m_node[j]);
   //    v[j] = node.get_vec3d(dofVIMX, dofVIMY, dofVIMZ) * m_alphaf + node.get_vec3d_prev(dofVIMX, dofVIMY, dofVIMZ) * (1.0-m_alphaf);
   // }
   // mat3d gradv;

   // strain rate
   mat3ds ev;

	// repeat for all integration points
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		double detJt = (m_update_dynamic ? invjact(el, Ji, n, m_alphaf) : invjact(el, Ji, n));

		detJt *= gw[n];

		const double* Gr = el.Gr(n);
		const double* Gs = el.Gs(n);
		const double* Gt = el.Gt(n);

      // // fluid viscous stress
      // gradv = gradient(el, v, n, m_alphaf);
      // ev = gradv.sym();

      // fluid viscous stress
      ev = pt.m_L.sym();

		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			double Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			double Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			double Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];        

			// calculate the force
			fe[3*i  ] -= detJt * ( 2.0 * fluid_dyn_viscosity * ev(0,0) * Gx 
                              + 2.0 * fluid_dyn_viscosity * ev(0,1) * Gy 
                              + 2.0 * fluid_dyn_viscosity * ev(0,2) * Gz);

			fe[3*i+1] -= detJt * ( 2.0 * fluid_dyn_viscosity * ev(1,0) * Gx 
                              + 2.0 * fluid_dyn_viscosity * ev(1,1) * Gy 
                              + 2.0 * fluid_dyn_viscosity * ev(1,2) * Gz);

			fe[3*i+2] -= detJt * ( 2.0 * fluid_dyn_viscosity * ev(2,0) * Gx 
                              + 2.0 * fluid_dyn_viscosity * ev(2,1) * Gy 
                              + 2.0 * fluid_dyn_viscosity * ev(2,2) * Gz);
		}
	}
}

//-----------------------------------------------------------------------------
// TANGENT
void MFEMImmersedElasticDomain::AssembleDomainStiffness(FEGlobalMatrix &K, double am, double av, double af, std::vector<MFEMBodyForce*> &bfs)
{
   // repeat over all solid elements
	int NE = Elements();
	
	// #pragma omp parallel for shared (NE)
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

      // get mfem submesh element id
      int submesh_id = -1;
      if (solid_submesh)
      {
         // serial
         submesh_id = solid_submesh->GetSubMeshElementFromParent(el.GetID()-1);
      }         
      else
      {
         // parallel
         submesh_id = el.GetID()-1;
         // submesh_id = el.GetID();
      }

      // mfem::out << "\n";
      // mfem::out << "Element id = " << submesh_id << std::endl;
      // auto curr_elem = dfes->GetFE(submesh_id);
      // auto trans = dfes->GetElementTransformation(submesh_id);
      // auto Nodes = curr_elem->GetNodes();
      // mfem::Vector curr_pos(3);
      // for(int j=0; j<Nodes.GetNPoints(); ++j)
      // {
      //    auto ip_node = Nodes.IntPoint(j);
      //    trans->SetIntPoint(&ip_node);
      //    trans->Transform(ip_node, curr_pos);
      //    mfem::out << std::setprecision(8) << "Node[" << j << "]=(" << curr_pos(0) << "," << curr_pos(1) << "," << curr_pos(2) << ")" << std::endl;
      // }
   

		if (el.isActive()) {

			// get the element's LM vector
			vector<int> lm;
         mfem::Array<int> vdofs;
         dfes->GetElementVDofs(submesh_id,vdofs);
         mfem::Array<int> vdofs_febio = TransformVDofs(vdofs,el.Nodes()); // reorder to match febio elem vector [xyz,xyz,xyz,...]
		   lm.assign(vdofs_febio.GetData(), vdofs_febio.GetData() + vdofs_febio.Size());

			// element stiffness matrix
			FEElementMatrix ke(el, lm);

			// create the element's stiffness matrix
			int ndof = 3 * el.Nodes();
			ke.resize(ndof, ndof);
			ke.zero();

         // dR/da
         ElementMassMatrix(el,ke,am);

         // dR/dv
         if (subtract_visc_fluid)
         {
            ElementArtificialFluidStiffness(el, ke, av);
         }

         // dR/du
         ElementSolidStiffness(iel, el, ke, af);
         ElementShapeDerivativeStiffness(el, ke, af, bfs);

         // // Open a file to write the matrix
         // std::ofstream outfile("matrix_ke.txt", std::ios::out);

         // mfem::Array<int> febio_map(ndof);
         // for (int i=0; i<el.Nodes(); ++i)
         // {
         //    febio_map[i] = 3*i;
         //    febio_map[i+el.Nodes()] = 3*i+1;
         //    febio_map[i+2*el.Nodes()] = 3*i+2;
         // }

         // // // febio_map.Print();

         // int precision_out = 6;
         // int width_out = precision_out + 6;

         // // Print the matrix ke to the file
         // for (int i = 0; i < ndof; ++i) 
         // {
         //    const int ii = febio_map[i];
         //    for (int j = 0; j < ndof; ++j) 
         //    {
         //       const int jj = febio_map[j];
         //       if (abs(ke(ii,jj)) >= 1e-12)
         //       {
         //          outfile << std::setw(width_out) << std::setprecision(precision_out) << ke(ii,jj);
         //       }
         //       else if (isnan(ke(ii,jj)))
         //       {
         //          outfile << std::setw(width_out) << std::setprecision(precision_out) << ke(ii,jj);
         //       }
         //       else
         //       {  
         //          outfile << std::setw(width_out) << std::setprecision(precision_out) << 0.0;
         //       }                  
         //       if (j < ndof - 1) outfile << " ";
         //    }
         //    outfile << std::endl;
         // }
         // outfile << std::endl; // Add an extra line between matrices

         // // Close the file
         // outfile.close();

         // // check if any entry is NAN
         // for (int i=0; i<ndof; ++i)
         // {
         //    const int ii = febio_map[i];
         //    for (int j=0; j<ndof; ++j)
         //    {
         //       const int jj = febio_map[j];
         //       if (isnan(ke(ii,jj)))
         //       {
         //          mfem::out << "-----------------------------------\n";
         //          mfem::out << "Element stiffness matrix has a NAN!" << std::endl;
         //          mfem::out << "Element id = " << submesh_id << std::endl;    
         //          mfem::out << "Entry " << ii << " " << jj << " = " << ke(ii,jj) << std::endl;
         //          mfem::out << "-----------------------------------\n";

         //       }
         //    }
         // }
         // double norm = ke.inf_norm();
         // if(isnan(norm))
         // {
         //    mfem::out << "Element stiffness matrix has a NAN!" << std::endl;
         // }



			// assemble element matrix into global stiffness matrix
         mfem::SparseMatrix& mfemK = *dynamic_cast<MFEMGlobalMatrix*>(&K)->GetMfemMatrix();
         // mfem::DenseMatrix mfem_elmat(vdofs_febio.Size());
         // mfem_elmat = 0.0;
         // for (int i = 0; i < vdofs_febio.Size(); ++i)
         // {
         //    for (int j = 0; j < vdofs_febio.Size(); ++j)
         //    {
         //       mfem_elmat(i,j) = ke(i,j);
         //    }
         // }
         mfem::DenseMatrix mfem_elmat(*ke, vdofs_febio.Size(), vdofs_febio.Size());
         mfemK.AddSubMatrix(vdofs_febio, vdofs_febio, mfem_elmat, 0);

			// K.Assemble(ke);
		}
	}
}

void MFEMImmersedElasticDomain::ElementMassMatrix(FESolidElement& el, matrix& ke, double scale)
{
   // Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();

   // gauss weights
   double* w = el.GaussWeights();

   for (int n = 0; n < nint; ++n)
   {
      FEMaterialPoint& mp = *el.GetMaterialPoint(n);
      FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
      double dens = m_pMat->Density(mp) - fluid_density*pt.m_J;
      double J0 = detJ0(el, n)*w[n];

      // set the shape function values
      double* H = el.H(n);

      // calculate stiffness component
      for (int i = 0; i < neln; ++i)
         for (int j = 0; j < neln; ++j)
         {
            double tmp = scale * dens * H[i] * H[j] * J0;
            ke[3*i  ][3*j  ] += tmp;
				ke[3*i+1][3*j+1] += tmp;
				ke[3*i+2][3*j+2] += tmp;
         }
   }
}

void MFEMImmersedElasticDomain::ElementArtificialFluidStiffness(FESolidElement& el, matrix& ke, double scale)
{
   // Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();

	// global derivatives of shape functions
	vec3d G[FEElement::MAX_NODES];

	// jacobian
	double detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian and shape function gradients
		detJt = ShapeGradient(el, n, G, m_alphaf)*gw[n]*scale;

		for (int na=0; na<neln; ++na)
         for (int nb=0; nb<neln; ++nb)
         {
            double kab_diag = (G[na]*(G[nb]))*detJt*fluid_dyn_viscosity;
            double fac = fluid_dyn_viscosity*detJt;

            ke[3*na  ][3*nb  ] -= kab_diag + fac * G[na].x * G[nb].x;
            ke[3*na  ][3*nb+1] -= fac * G[na].y * G[nb].x;
            ke[3*na  ][3*nb+2] -= fac * G[na].z * G[nb].x;

            ke[3*na+1][3*nb  ] -= fac * G[na].x * G[nb].y;
            ke[3*na+1][3*nb+1] -= kab_diag + fac * G[na].y * G[nb].y;
            ke[3*na+1][3*nb+2] -= fac * G[na].z * G[nb].y;

            ke[3*na+2][3*nb  ] -= fac * G[na].x * G[nb].z;
            ke[3*na+2][3*nb+1] -= fac * G[na].y * G[nb].z;
            ke[3*na+2][3*nb+2] -= kab_diag + fac * G[na].z * G[nb].z;
         }
	}
}

void MFEMImmersedElasticDomain::ElementSolidStiffness(int iel, FESolidElement& el, matrix& ke, double scale)
{
   ElementGeometricalStiffness(iel, el, ke, scale);
   ElementMaterialStiffness(iel, el, ke, scale);
}

//! calculates element's geometrical stiffness component for integration point n
void MFEMImmersedElasticDomain::ElementGeometricalStiffness(int iel, FESolidElement &el, matrix &ke, double scale)
{
	// spatial derivatives of shape functions
	vec3d G[FEElement::MAX_NODES];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate geometrical element stiffness matrix
	int neln = el.Nodes();
	int nint = el.GaussPoints();
	for (int n = 0; n<nint; ++n)
	{
		// calculate shape function gradients and jacobian
      double w = ShapeGradient(el, n, G, m_alphaf)*gw[n]*scale;

		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// element's Cauchy-stress tensor at gauss point n
		mat3ds& s = pt.m_s;

		for (int i = 0; i<neln; ++i)
			for (int j = 0; j<neln; ++j)
			{
				double kab = (G[i]*(s * G[j]))*w;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

void MFEMImmersedElasticDomain::ElementMaterialStiffness(int iel, FESolidElement &el, matrix &ke, double scale)
{
	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();

	// global derivatives of shape functions
	vec3d G[FEElement::MAX_NODES];

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	// jacobian
	double detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian and shape function gradients
		detJt = ShapeGradient(el, n, G, m_alphaf)*gw[n]*scale;

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		// get the 'D' matrix
//		tens4ds C = m_pMat->Tangent(mp);
        tens4dmm C = (m_secant_tangent ? m_pMat->SecantTangent(mp) : m_pMat->SolidTangent(mp));
		C.extract(D);

		for (int i=0, i3=0; i<neln; ++i, i3 += 3)
		{
			Gxi = G[i].x;
			Gyi = G[i].y;
			Gzi = G[i].z;

			for (int j=0, j3 = 0; j<neln; ++j, j3 += 3)
			{
				Gxj = G[j].x;
				Gyj = G[j].y;
				Gzj = G[j].z;

				// calculate D*BL matrices
				DBL[0][0] = (D[0][0]*Gxj+D[0][3]*Gyj+D[0][5]*Gzj);
				DBL[0][1] = (D[0][1]*Gyj+D[0][3]*Gxj+D[0][4]*Gzj);
				DBL[0][2] = (D[0][2]*Gzj+D[0][4]*Gyj+D[0][5]*Gxj);

				DBL[1][0] = (D[1][0]*Gxj+D[1][3]*Gyj+D[1][5]*Gzj);
				DBL[1][1] = (D[1][1]*Gyj+D[1][3]*Gxj+D[1][4]*Gzj);
				DBL[1][2] = (D[1][2]*Gzj+D[1][4]*Gyj+D[1][5]*Gxj);

				DBL[2][0] = (D[2][0]*Gxj+D[2][3]*Gyj+D[2][5]*Gzj);
				DBL[2][1] = (D[2][1]*Gyj+D[2][3]*Gxj+D[2][4]*Gzj);
				DBL[2][2] = (D[2][2]*Gzj+D[2][4]*Gyj+D[2][5]*Gxj);

				DBL[3][0] = (D[3][0]*Gxj+D[3][3]*Gyj+D[3][5]*Gzj);
				DBL[3][1] = (D[3][1]*Gyj+D[3][3]*Gxj+D[3][4]*Gzj);
				DBL[3][2] = (D[3][2]*Gzj+D[3][4]*Gyj+D[3][5]*Gxj);

				DBL[4][0] = (D[4][0]*Gxj+D[4][3]*Gyj+D[4][5]*Gzj);
				DBL[4][1] = (D[4][1]*Gyj+D[4][3]*Gxj+D[4][4]*Gzj);
				DBL[4][2] = (D[4][2]*Gzj+D[4][4]*Gyj+D[4][5]*Gxj);

				DBL[5][0] = (D[5][0]*Gxj+D[5][3]*Gyj+D[5][5]*Gzj);
				DBL[5][1] = (D[5][1]*Gyj+D[5][3]*Gxj+D[5][4]*Gzj);
				DBL[5][2] = (D[5][2]*Gzj+D[5][4]*Gyj+D[5][5]*Gxj);

				ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*detJt;
				ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*detJt;
				ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*detJt;

				ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*detJt;
				ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*detJt;
				ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*detJt;

				ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*detJt;
				ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*detJt;
				ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*detJt;
			}
		}
	}
}

void MFEMImmersedElasticDomain::ElementShapeDerivativeStiffness(FESolidElement& el, matrix& ke, double scale, std::vector<MFEMBodyForce*>& bfs)
{
   // Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();

	// global derivatives of shape functions in current configuration
	vec3d Gx[FEElement::MAX_NODES];

   // global derivatives of shape functions in reference configuration
	vec3d GX[FEElement::MAX_NODES];

   // shape functions
   double *H;

   // strain rate
   mat3ds ev;

   // deformation gradient
   mat3d F, FinvT;

   // inertial force
   vec3d inertial_force, rho_if;

   // jacobian of wrt ref config 
   double detJ0;

	// jacobian
	double detJt;
	
	// weights at gauss points
	const double *gw = el.GaussWeights();

	// calculate element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
      FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate jacobian and shape function gradients
		detJt = ShapeGradient(el, n, Gx, m_alphaf)*gw[n]*scale;
      detJ0 = ShapeGradient0(el, n, GX);
      H = el.H(n);

      // fluid viscous stress
      ev = pt.m_L.sym();

      // deformation gradient inv transpose
      F = pt.m_F;
      FinvT = F.transinv();

      // inertial force
      // acceleration
      inertial_force(0) = pt.m_a(0);
      inertial_force(1) = pt.m_a(1);
      inertial_force(2) = pt.m_a(2);
      // body forces
      for (auto &bf : bfs)
      {
         // get the force
         inertial_force(0) -= bf->m_bx(mp);
         inertial_force(1) -= bf->m_by(mp);
         inertial_force(2) -= bf->m_bz(mp);
      }
      rho_if = inertial_force*fluid_density; // rhof * (a - f)

		for (int na=0; na<neln; ++na)
         for (int nb=0; nb<neln; ++nb)
         {
            double fac1 = detJt * H[na];
            double fac2 = subtract_visc_fluid ? detJt * 2.0 * fluid_dyn_viscosity : 0.0;

            double FinvTGX_x = FinvT(0,0)*GX[nb].x + FinvT(0,1)*GX[nb].y + FinvT(0,2)*GX[nb].z;
            double FinvTGX_y = FinvT(1,0)*GX[nb].x + FinvT(1,1)*GX[nb].y + FinvT(1,2)*GX[nb].z;
            double FinvTGX_z = FinvT(2,0)*GX[nb].x + FinvT(2,1)*GX[nb].y + FinvT(2,2)*GX[nb].z;

            double ev_Gx_x = ev(0,0) * Gx[na].x + ev(0,1) * Gx[na].y + ev(0,2) * Gx[na].z;
            double ev_Gx_y = ev(1,0) * Gx[na].x + ev(1,1) * Gx[na].y + ev(1,2) * Gx[na].z;
            double ev_Gx_z = ev(2,0) * Gx[na].x + ev(2,1) * Gx[na].y + ev(2,2) * Gx[na].z;

            ke[3*na  ][3*nb  ] -= fac1 * rho_if.x * FinvTGX_x + fac2 * ev_Gx_x * FinvTGX_x;
            ke[3*na  ][3*nb+1] -= fac1 * rho_if.x * FinvTGX_y + fac2 * ev_Gx_x * FinvTGX_y;
            ke[3*na  ][3*nb+2] -= fac1 * rho_if.x * FinvTGX_z + fac2 * ev_Gx_x * FinvTGX_z;

            ke[3*na+1][3*nb  ] -= fac1 * rho_if.y * FinvTGX_x + fac2 * ev_Gx_y * FinvTGX_x;
            ke[3*na+1][3*nb+1] -= fac1 * rho_if.y * FinvTGX_y + fac2 * ev_Gx_y * FinvTGX_y;
            ke[3*na+1][3*nb+2] -= fac1 * rho_if.y * FinvTGX_z + fac2 * ev_Gx_y * FinvTGX_z;

            ke[3*na+2][3*nb  ] -= fac1 * rho_if.z * FinvTGX_x + fac2 * ev_Gx_z * FinvTGX_x;
            ke[3*na+2][3*nb+1] -= fac1 * rho_if.z * FinvTGX_y + fac2 * ev_Gx_z * FinvTGX_y;
            ke[3*na+2][3*nb+2] -= fac1 * rho_if.z * FinvTGX_z + fac2 * ev_Gx_z * FinvTGX_z;
         }
	}
}

//-----------------------------------------------------------------------------
void MFEMImmersedElasticDomain::BuildMatrixProfile(FEGlobalMatrix& M)
{
	vector<int> elm;
   mfem::Array<int> vdofs;
	const int NE = Elements();
	for (int j = 0; j<NE; ++j)
	{
      FESolidElement& el = m_Elem[j];
      int submesh_id = solid_submesh ? solid_submesh->GetSubMeshElementFromParent(el.GetID()-1) : el.GetID() - 1;
      dfes->GetElementVDofs(submesh_id,vdofs);
      elm.assign(vdofs.GetData(), vdofs.GetData() + vdofs.Size()); // check this
		M.build_add(elm);
	}
}