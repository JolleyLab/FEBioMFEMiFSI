#pragma once
#include <ImmersedFSI-Serial/immersed_solid.hpp>
#include <FEBioMech/FEElasticSolidDomain.h>
#include "MFEMLoad.h"
#include "MFEMTools.h"
#include <FEBioMech/FEElasticMaterial.h>

using namespace immersed_solid_serial;

class MFEMImmersedElasticDomain : public FESolidDomain
{
public:
   // constructor
   MFEMImmersedElasticDomain(FEModel* pfem);

   // destructor
   virtual ~MFEMImmersedElasticDomain() noexcept;   

   // set FESpace
   void SetFESpace(mfem::FiniteElementSpace *fes)
   {
      dfes = fes;
   };

   // get the total dof list
   const FEDofList& GetDOFList() const override
   {
      return m_dof;
   };

   // set fluid properties
   void SetFluidProperties(double rho, double mu);

   // set additional time stepping data
   void SetAdditionalData(bool subtract_visc_fluid_, mfem::SubMesh* solid_submesh_, mfem::FiniteElementSpace *gfes_);

   void SetMaterial(FEMaterial* pmat) override
	{
      if (dynamic_cast<FEElasticMaterial*>(pmat)) 
      {
         m_pMat = dynamic_cast<FEElasticMaterial*>(pmat);
         FESolidDomain::SetMaterial(pmat);
      }
      else assert(false);
	}

	FEMaterial* GetMaterial() override { return m_pMat; }

   mfem::Array<int> &TransformVDofs(mfem::Array<int> &vdofs, int nd);

public:
   //! initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

   // update stresses
	void Update(const FETimeInfo& tp) override;

	// update the element stress
	virtual void UpdateElementStress(int iel, const FETimeInfo& tp);
   
   //! activate
	void Activate() override;

   // unpack lm data
   void UnpackLM(FEElement& el, vector<int>& lm) override;

public:
	// --- S T I F F N E S S ---
   void AssembleDomainStiffness(FEGlobalMatrix &K, double am, double av, double af, std::vector<MFEMBodyForce*> &bfs);

   // mass matrix
   void ElementMassMatrix(FESolidElement& el, matrix& ke, double scale);

   // artificial fluid viscous stress
   void ElementArtificialFluidStiffness(FESolidElement& el, matrix& ke, double scale);

   // solid stiffness
   virtual void ElementSolidStiffness(int iel, FESolidElement& el, matrix& ke, double scale);

	//! geometrical stiffness (i.e. initial stress)
	virtual void ElementGeometricalStiffness(int iel, FESolidElement& el, matrix& ke);

	//! material stiffness component
	virtual void ElementMaterialStiffness(int iel, FESolidElement& el, matrix& ke);

   // shape derivative stiffness
   void ElementShapeDerivativeStiffness(FESolidElement& el, matrix& ke, double scale, std::vector<MFEMBodyForce*>& bfs);
   
	// --- R E S I D U A L ---
   void AssembleDomainResidual(MFEMGlobalVector& R, std::vector<MFEMBodyForce*> &bfs);

   //! Calculates the internal stress vector for solid elements
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

   //! Calculates the inertial force vector for solid elements
   void ElementInertialForce(FESolidElement& el, vector<double>& fe, std::vector<MFEMBodyForce*> &bfs);

   //! Calculates the artificial fluid force vector
   void ElementArtificialFluidForce(FESolidElement& el, vector<double>& fe);

public:
   //! build the matrix profile only for the solid domain in the mesh
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

protected:
   FEDofList	m_dof;		// total dof list

   // fluid properties (for now constants)
   double fluid_density; // density
   double fluid_dyn_viscosity; // dynamic viscosity
   bool subtract_visc_fluid = false; // flag for subtracting artificial fluid viscous stress

   // additional time stepping factors
   double dt;
   double m_gamma;
   double m_alphaf;
   double m_alpham;
   double m_beta;

   // fespace defined by MFEM (equivalent numbering to FEBio)
   mfem::FiniteElementSpace *dfes = nullptr; // local fespace
   mfem::FiniteElementSpace *gfes = nullptr; // global fespace

   // solid submesh from MFEM
   mfem::SubMesh *solid_submesh = nullptr;

   // additional data
   bool m_update_dynamic = true;
   bool m_secant_tangent = false;
   bool m_secant_stress = false;
   mfem::Vector mfem_fe;
   mfem::Array<int> febio_vdofs;

   // material
   FESolidMaterial*	m_pMat;
};


