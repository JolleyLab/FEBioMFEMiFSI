#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEBodyLoad.h>
#include <FECore/FEModelParam.h>
#include <FECore/FENodalLoad.h>

class MFEMSolver;
namespace mfem
{
   class Mesh;
}

/////// Base classes ///////

class MFEMSurfaceLoad : public FESurfaceLoad
{
public:
   MFEMSurfaceLoad(FEModel* fem) : FESurfaceLoad(fem), m_backflow_stab(false) 
   {
      m_beta = 0.2;
   };

   // Get the mesh
   mfem::Mesh* GetMesh(MFEMSolver* solver);

	virtual void Apply(MFEMSolver* solver, int attr);

   virtual void CopyFrom(FEModelLoad* pml);

   void SetBoundaryID(int id) { m_bdr_attr = id; };

   double GetBoundaryID() { return m_bdr_attr; };

private:
   int m_bdr_attr; // bdr attribute ID for the boundary condition
   bool m_backflow_stab; // flag for backflow stabilization
	FEParamDouble m_beta;

   DECLARE_FECORE_CLASS();
};

class MFEMBodyLoad : public FEBodyLoad
{
public:
	MFEMBodyLoad(FEModel* fem) : FEBodyLoad(fem) {}

   // Get the mesh
   mfem::Mesh* GetMesh(MFEMSolver* solver);

	virtual void Apply(MFEMSolver* solver, int attr) = 0;

   void SetDomainIDList(std::vector<int> &list) { m_dom_list = list; };

   std::vector<int> &GetDomainIDList() { return m_dom_list; };

   virtual void CopyFrom(FEModelLoad* pml);

private:
   std::vector<int> m_dom_list; // domain IDs for the body load
};

class MFEMNodalLoad : public FENodalLoad
{
public:
	MFEMNodalLoad(FEModel* fem) : FENodalLoad(fem) {}

   // override pure virtual functions
   void GetNodalValues(int inode, std::vector<double>& vals) override {};
   bool SetDofList(FEDofList& doflist) override 
   { 
     return doflist.AddDof(m_dof);
   };

	virtual void Apply(MFEMSolver* solver) = 0;

   // virtual void CopyFrom(FEModelLoad* pml);

protected:
   int m_dof;
};

/////// Derived classes for each bc ///////

class MFEMNeumannBC : public MFEMSurfaceLoad
{
public:
	MFEMNeumannBC(FEModel* fem) : MFEMSurfaceLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEModelLoad* pml) override;

private:
	FEParamDouble m_tx;
	FEParamDouble m_ty;
	FEParamDouble m_tz;

	DECLARE_FECORE_CLASS();
};

class MFEMTractionFreeBC : public MFEMSurfaceLoad
{
public:
	MFEMTractionFreeBC(FEModel* fem) : MFEMSurfaceLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEModelLoad* pml) override;

private:

	DECLARE_FECORE_CLASS();
};

class MFEMPressureBC : public MFEMSurfaceLoad
{
public:
	MFEMPressureBC(FEModel* fem) : MFEMSurfaceLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEModelLoad* pml) override;

private:
   // applies p*nor to the surface
	FEParamDouble m_p;

	DECLARE_FECORE_CLASS();
};

// class MFEMBackflowStabilizationBC : public MFEMSurfaceLoad
// {
// public:
// 	MFEMBackflowStabilizationBC(FEModel* fem);

// 	void Apply(MFEMSolver* solver, int attr) override;

//    void CopyFrom(FEModelLoad* pml) override;

// private:
//    // assumes constant beta
// 	FEParamDouble m_beta;

// 	DECLARE_FECORE_CLASS();
// };

class MFEMResistanceBC : public MFEMSurfaceLoad
{
public:
	MFEMResistanceBC(FEModel* fem) : MFEMSurfaceLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEModelLoad* pml) override;

private:
   // assumes constant R and p0
	FEParamDouble m_R; // proximal resistance
	FEParamDouble m_p0; // initial pressure

	DECLARE_FECORE_CLASS();
};

class MFEMWindkesselBC : public MFEMSurfaceLoad
{
public:
	MFEMWindkesselBC(FEModel* fem) : MFEMSurfaceLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEModelLoad* pml) override;

private:
   // assumes constant R, Rd, C, p0(spatial), pd
	FEParamDouble m_R; // proximal resistance 
	FEParamDouble m_Rd; // distal resistance
   FEParamDouble m_C; // vascular compliance
   FEParamDouble m_p0; // initial pressure // TODO make this spatially dependent thru a coeff
   FEParamDouble m_pd; // distal pressure
   FEParamDouble m_Q0; // initial flow rate

	DECLARE_FECORE_CLASS();
};

class MFEMBodyForce : public MFEMBodyLoad
{
public:
	MFEMBodyForce(FEModel* fem) : MFEMBodyLoad(fem) {}

	void Apply(MFEMSolver* solver, int attr) override;

   void CopyFrom(FEModelLoad* pml) override;

public:
   FEParamDouble m_bx;
	FEParamDouble m_by;
	FEParamDouble m_bz;

	DECLARE_FECORE_CLASS();
};

class MFEMSetNodalPressure : public MFEMNodalLoad
{
public:
	MFEMSetNodalPressure(FEModel* fem) : MFEMNodalLoad(fem) {}

	void Apply(MFEMSolver* solver) override;

   // void CopyFrom(FEModelLoad* pml) override;
   
private:
   FEParamDouble m_p;

	DECLARE_FECORE_CLASS();
};