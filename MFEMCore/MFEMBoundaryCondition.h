#pragma once
#include <FECore/FESurfaceBC.h>
#include <FECore/FEModelParam.h>

class MFEMSolver;
namespace mfem
{
   class Mesh;
}

// Base class for boundary conditions that can be used in the MFEMSolver
class MFEMBoundaryCondition : public FESurfaceBC
{
public:
	MFEMBoundaryCondition(FEModel* fem) : FESurfaceBC(fem) {};

	virtual void CopyFrom(FEBoundaryCondition* pbc) override;

   // Get the mesh
   mfem::Mesh* GetMesh(MFEMSolver* solver);

	// This function must be overridden by derived classes and should be used
	// to apply the boundary condition to the MFEMSolver
	virtual void Apply(MFEMSolver* solver, int attr) = 0;

   void SetBoundaryID(int id) { m_bdr_attr = id; };
   int GetBoundaryID() { return m_bdr_attr; };

   void SetSurfaceName(const std::string &name) { m_surfname = name; };
   const std::string& GetSurfaceName() const { return m_surfname; };

   bool Init() override;

private:
   int m_bdr_attr; // bdr attribute ID for the boundary condition
   std::string m_surfname;
};

// Dirichlet velocity boundary condition
class MFEMVelDirichletBC : public MFEMBoundaryCondition
{
public:
	MFEMVelDirichletBC(FEModel* fem) : MFEMBoundaryCondition(fem) {};

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEBoundaryCondition* pbc) override;

private:
	FEParamDouble m_vx;
	FEParamDouble m_vy;
	FEParamDouble m_vz;

   bool m_zero_bdr = false; // if true, the boundary nodes of this bc are set to zero velocity

	DECLARE_FECORE_CLASS();
};

// No slip bc
class MFEMNoSlipBC : public MFEMBoundaryCondition
{
public:
	MFEMNoSlipBC(FEModel* fem) : MFEMBoundaryCondition(fem) {};

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEBoundaryCondition* pbc) override;

private:

	DECLARE_FECORE_CLASS();
};

class MFEMPrescribedVelocityComponentBC : public MFEMBoundaryCondition
{
public:
	MFEMPrescribedVelocityComponentBC(FEModel* fem) : MFEMBoundaryCondition(fem) {};

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEBoundaryCondition* pbc) override;

private:
   FEParamDouble m_vval;
   int m_dof;
   bool m_zero_bdr = false; // if true, the boundary nodes of this bc are set to zero velocity

	DECLARE_FECORE_CLASS();
};

// Set periodic bcs (on both velocity and pressure)
class MFEMPeriodicBC : public MFEMBoundaryCondition
{
public:
	MFEMPeriodicBC(FEModel* fem) : MFEMBoundaryCondition(fem) {};

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEBoundaryCondition* pbc) override;

private:

	DECLARE_FECORE_CLASS();
};

class MFEMMappedInflowBC : public MFEMBoundaryCondition
{
public:
	MFEMMappedInflowBC(FEModel* fem) : MFEMBoundaryCondition(fem) { m_profile = 0;};

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEBoundaryCondition* pbc) override;

private:
   std::string m_flowfile;
   bool m_periodic;
   int m_profile; // default is 0 (parabolic)

	DECLARE_FECORE_CLASS();
};