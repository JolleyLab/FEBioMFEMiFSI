#include "MFEMGlobalTypes.h"
#include <FECore/FEModel.h>
#include "MFEMImmersedElasticDomain.h"
#include <mfem.hpp>


MFEMGlobalVector::MFEMGlobalVector(FEModel* fem, int size, std::vector<double>& R, std::vector<double>& Fr) : FEGlobalVector(*fem, R, Fr)
{
   mfem_R = new mfem::Vector(R.data(),size);
}

MFEMGlobalVector::~MFEMGlobalVector() noexcept
{
   delete mfem_R;
}

//-----------------------------------------------------------------------------
void MFEMGlobalVector::Assemble(vector<int>& en, vector<int>& elm, vector<double>& fe, bool bdom)
{
	vector<double>& R = m_R;

	// assemble the element residual into the global residual
	int ndof = (int)fe.size();
	for (int i=0; i<ndof; ++i)
	{
		int I = elm[i];
		if ( I >= 0) {
// #pragma omp atomic
			R[I] += fe[i];
		}
	}
}

//-----------------------------------------------------------------------------
//! \todo This function does not add to m_Fr. Is this a problem?
void MFEMGlobalVector::Assemble(vector<int>& lm, vector<double>& fe)
{
	vector<double>& R = m_R;
	const int n = (int) lm.size();
	for (int i=0; i<n; ++i)
	{
		int nid = lm[i];
		if (nid >= 0) {
// #pragma omp atomic
			R[nid] += fe[i];
		}
	}
}

//-----------------------------------------------------------------------------
//! assemble a nodel value
void MFEMGlobalVector::Assemble(int nodeId, int dof, double f)
{
	// get the equation number
	FENode& node = m_fem.GetMesh().Node(nodeId);
	int n = node.m_ID[dof];

	// assemble into global vector
	if (n >= 0) {
// #pragma omp atomic
		m_R[n] += f;
	}
}

void MFEMGlobalVector::Assemble(mfem::Array<int> &vdofs, std::vector<double>& fe)
{
   mfem::Vector mfem_fe(fe.data(),fe.size());
   mfem_R->AddElementVector(vdofs, mfem_fe);
}

void MFEMGlobalVector::Assemble(mfem::Array<int> &vdofs, mfem::Vector &fe)
{
   mfem_R->AddElementVector(vdofs, fe);
}

MFEMGlobalMatrix::MFEMGlobalMatrix(CRSSparseMatrix* pK, bool del) : FEGlobalMatrix(pK, del)
{

}

void MFEMGlobalMatrix::CreateProfile(FEModel* pfem, int neq, bool reshape)
{
   build_begin(neq);
   {
      // setup static profile
      if(reshape)
      {
         m_MPs.Clear();

         // build the profile based on each domain
         (*this).BuildMatrixProfile(pfem, true);

         build_flush();
			m_MPs = *m_pMP;
      }
      else
      {
         // copy the old static profile
			*m_pMP = m_MPs;
      }

      // build the dynamic profile
      (*this).BuildMatrixProfile(pfem, false);      
   }
   build_end();

   // wrap FEBio matrix into MFEM SparseMatrix
   if(mfem_K) delete mfem_K;
   int *I = m_pA->Pointers();
   int *J = m_pA->Indices();
   double *data = m_pA->Values();
   mfem_K = new mfem::SparseMatrix(I,J,data,neq,neq,false,false,true);
}

void MFEMGlobalMatrix::BuildMatrixProfile(FEModel* pfem, bool breset)
{
   FEModel& fem = *pfem;
	FEMesh& mesh = fem.GetMesh();

	// when reset is true we build the entire matrix profile
	// (otherwise we only build the "dynamic" profile)
	if (breset)
	{
		// Add all elements in the immersed solid domain to the profile      
		for (int nd = 0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& d = mesh.Domain(nd);
         if (dynamic_cast<MFEMImmersedElasticDomain*>(&d))
         {
            d.BuildMatrixProfile(*this);
         }			
		}

		// // linear constraints
		// FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
		// LCM.BuildMatrixProfile(*this);
	}
	else
	{
      
      // TODO: contact etc


		// // Do the "dynamic" profile. That is the part of the profile that always changes
		// // This is mostly contact
		// // do the nonlinear constraints
		// int M = fem.NonlinearConstraints();
		// for (int m = 0; m<M; ++m)
		// {
		// 	FENLConstraint* pnlc = fem.NonlinearConstraint(m);
		// 	if (pnlc->IsActive()) pnlc->BuildMatrixProfile(*this);
		// }

		// // All following "elements" are nonstatic. That is, they can change
		// // connectivity between calls to this function. All of these elements
		// // are related to contact analysis (at this point).
		// if (fem.SurfacePairConstraints() > 0)
		// {
		// 	// Add all contact interface elements
		// 	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
		// 	{
		// 		FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
		// 		if (pci->IsActive()) pci->BuildMatrixProfile(*this);
		// 	}
		// }
	}
}