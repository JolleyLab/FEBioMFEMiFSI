#pragma once
#define IFS_VERSION 1.0

#include "mfem.hpp"

#include "../IncompressibleFluid/ifs_solver_enum.hpp"
#include "../IncompressibleFluid/ifs_bc_containers.hpp"

#include "SpatialOperators/ResidualOpFluid.hpp"
#include "SpatialOperators/TractionIntegrator.hpp"
#include "SpatialOperators/BodyForceIntegrator.hpp"
#include "SpatialOperators/IncompressibleFluidIntegrator.hpp"

#include "CardiovascularBCs/BackflowStabIntegrator.hpp"
#include "CardiovascularBCs/CVBC_Op.hpp"
#include "CardiovascularBCs/mapped_inflow_bc.hpp"

#include "../ImmersedFSI/ImmersedSolid.hpp"

#include "Solvers/BlockOrthoSolver.hpp"
#include "Solvers/FluidNewtonSolver.hpp"
#include "Solvers/Preconditioners/BlockPreconditioner.hpp"
#include "Solvers/Preconditioners/SIMPLE_Prec.hpp"
#include "Solvers/Preconditioners/SCR_Prec.hpp"
#include "Solvers/Preconditioners/aSCR_Prec.hpp"
#include "Solvers/Preconditioners/CVBC_AinvPC.hpp"

#include "../Utilities/extended_solvers.hpp"
#include "../Utilities/MinItMonitor.hpp"
#include "../Utilities/TimerRegistry.hpp"

namespace ifs_solver
{
// Class for Incompressible IncompressibleFluid Solver based on CG Method
class IncompressibleFluidSolver
{
public:

   ////////////////////
   /// CONSTRUCTORS ///
   ////////////////////

   /// Initialize data structures, set FE space order and kinematic viscosity.
   #ifdef MFEM_USE_MPI
   IncompressibleFluidSolver(mfem::ParMesh *mesh, int orderv, int orderp, double kin_vis, int ordering_type_ = mfem::Ordering::byNODES);
   // constructor for dimensional form of equations
   IncompressibleFluidSolver(mfem::ParMesh *mesh, int orderv, int orderp, double kin_vis, double fluid_density_, int ordering_type_ = mfem::Ordering::byNODES) : IncompressibleFluidSolver(mesh, orderv, orderp, kin_vis, ordering_type_)
   {
      fluid_density = fluid_density_;
      rhof.constant = fluid_density;
      // dynamic_viscosity.constant = kin_vis * fluid_density;
      mu_f.constant = kin_vis * fluid_density;
      dimensional_flag = ifs_solver::dim_flag::Dimensional;
   };
   #else
   IncompressibleFluidSolver(mfem::Mesh *mesh, int orderv, int orderp, double kin_vis, int ordering_type_ = mfem::Ordering::byNODES);
   // constructor for dimensional form of equations
   IncompressibleFluidSolver(mfem::Mesh *mesh, int orderv, int orderp, double kin_vis, double fluid_density_, int ordering_type_ = mfem::Ordering::byNODES) : IncompressibleFluidSolver(mesh, orderv, orderp, kin_vis, ordering_type_)
   {
      fluid_density = fluid_density_;
      rhof.constant = fluid_density;
      // dynamic_viscosity.constant = kin_vis * fluid_density;
      mu_f.constant = kin_vis * fluid_density;
      dimensional_flag = ifs_solver::dim_flag::Dimensional;
   };
   #endif
   

   /////////////////////////
   /// DRIVING FUNCTIONS ///
   /////////////////////////

   /// Initialize forms, solvers and preconditioners.
   void Setup(double dt);

   // Perform a time step (time = t_n+1) so the variables stored in
   // un_gf and pn_gf correspond to the current time after this function
   // is called
   void Step(double &time, double dt, int cur_step, bool provisional);
   bool Step(double &time, double dt, int cur_step);

   /// Remove mean from a Vector.
   /**
    * Modify the Vector @a v by subtracting its mean using
    * \f$v = v - \frac{\sum_i^N v_i}{N} \f$
    */
   void Orthogonalize(mfem::Vector &v);

   /// Remove the mean from a GridFunction.
   /**
    * Modify the ParGridFunction @a v by subtracting its mean using
    * \f$ v = v - \int_\Omega \frac{v}{vol(\Omega)} dx \f$.
    */
   void MeanZero(mfem::GridFunction &v);

   /// Rotate entries in the time step and solution history arrays.
   void UpdateTimestepHistory(const double &time,double &dt);
   
   /// Print timing summary of the solving routine.
   virtual void PrintTimingData();
   void ConvertTimeToFormattedTime(double time_in_seconds, char* time_str, size_t size);
   mfem::StopWatch sw_io, sw_step;

   void UpdateParameters(double time, double dt, int current_step_);
   void Predictors(double dt);
   void SetEssentialBCs(const double &time);
   void UpdateVariables(double &time, double &dt);
   void UpdateWindkesselParameters();
   void GetNewtonSolverConvergenceInfo(int &newton_its, double &newton_res, double &newton_relres)
   {
      newton_its = newton_solver->GetNumIterations();
      newton_res = newton_solver->GetFinalNorm();
      newton_relres = newton_solver->GetFinalRelNorm();
   };

   //////////////////////////
   /// GET DATA FUNCTIONS ///
   //////////////////////////

   /// Return a pointer to the current velocity ParGridFunction.
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction *GetCurrentVelocity() { return &vn_gf; }
   #else
   mfem::GridFunction *GetCurrentVelocity() { return &vn_gf; }
   #endif

   /// Return a pointer to the current pressure ParGridFunction.
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction *GetCurrentPressure() { return &pn_gf; }
   #else
   mfem::GridFunction *GetCurrentPressure() { return &pn_gf; }
   #endif

   /// Return a pointer to the immersed solid displacement
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction *GetImmersedSolidDisplacement() { return im_solid->GetCurrentDisplacement(); }
   #else
   mfem::GridFunction *GetImmersedSolidDisplacement() { return im_solid->GetCurrentDisplacement(); }
   #endif

   /// Return a pointer to the immersed solid velocity
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction *GetImmersedSolidVelocity() { return im_solid->GetCurrentVelocity(); }
   #else
   mfem::GridFunction *GetImmersedSolidVelocity() { return im_solid->GetCurrentVelocity(); }
   #endif

   /// Return a pointer to the immersed solid velocity
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction *GetImmersedSolidAcceleration() { return im_solid->GetCurrentAcceleration(); }
   #else
   mfem::GridFunction *GetImmersedSolidAcceleration() { return im_solid->GetCurrentAcceleration(); }
   #endif

   /// Return a pointer to the immersed solid
   immersed_solid::ImmersedSolid *GetImmersedSolid() { return im_solid; }

   /// Get the indicator function for the immersed solid boundary
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction *GetImmersedSolidIndicator() { return &s_gf; }
   #else
   mfem::GridFunction *GetImmersedSolidIndicator() { return &s_gf; }
   #endif

   // Get the mesh
   #ifdef MFEM_USE_MPI
   mfem::ParMesh *GetMesh() { return pmesh; }
   #else
   mfem::Mesh *GetMesh() { return pmesh; }
   #endif


   ///////////////////////////
   /// BOUNDARY CONDITIONS ///
   ///////////////////////////

   //////////////// Velocity BCs ////////////////
   /// Add a Dirichlet boundary condition to the velocity field in cartesian components
   void AddVelDirichletBC(mfem::VectorCoefficient *coeff, mfem::Array<int> &attr, bool zero_bdr = false);
   void AddVelDirichletBC(VecFuncT *f, mfem::Array<int> &attr, bool zero_bdr = false);
   // components are 1,2,3 = x,y,z
   void AddVelDirichletBC(mfem::Coefficient *coeff, mfem::Array<int> &attr,int component_, bool zero_bdr = false);
   void AddVelDirichletBC(ScalarFuncT *f, mfem::Array<int> &attr,int component_, bool zero_bdr = false);

   // Add no slip velocity boundary condition
   void AddNoSlipBC(mfem::Array<int> &attr);
   static void no_slip(const mfem::Vector &x, double t, mfem::Vector &u)
   {
      u = 0.0;
   };

   //////////////// Periodic BCs ////////////////
   // Add periodic boundary condition to both velocity and pressure
   void AddPeriodicBC(mfem::Array<int> &attr);

   //////////////// Traction BCs ////////////////
   /// Add a neumann boundary condition to the velocity field.
   void AddVelNeumannBC(mfem::VectorCoefficient *coeff, mfem::Array<int> &attr); // here user specifies gradu * n entirely
   void AddVelNeumannBC(VecFuncT *f, mfem::Array<int> &attr);

   // overload to pass grad_v directly, so here bc is specified as gradu (user determined) * n
   void AddVelNeumannBC(mfem::MatrixCoefficient *coeff, mfem::Array<int> &attr);
   void AddVelNeumannBC(MatFuncT *f, mfem::Array<int> &attr);

   // overload for scalar pressure coefficient
   void AddVelNeumannBC(mfem::Coefficient *coeff, mfem::Array<int> &attr, bool pressure = false); // here user specifies gradu * n entirely
   void AddVelNeumannBC(ScalarFuncT *f, mfem::Array<int> &attr, bool pressure = false);

   // Add a pressure bc as normal traction, this is essentially AddVelNeumannBC with zero_mean = false
   /*
      Assumes specified pressure is in same direction as normal on surface. To get a compressive force
      set scalar to -p! Can also use wrapper AddPressureBC to handle this directly
   */
   void AddNormalTractionBC(mfem::Coefficient *coeff, mfem::Array<int> &attr, bool pressure = false)
   {
      MFEM_ASSERT(viscous_term == ifs_solver::ViscousTermType::StressForm, "Normal traction BC only supported for stress form of viscous term");
      zero_mean = false;
      this->AddVelNeumannBC(coeff,attr,pressure);
   }
   void AddNormalTractionBC(ScalarFuncT *f, mfem::Array<int> &attr, bool pressure = false)
   {
      MFEM_ASSERT(viscous_term == ifs_solver::ViscousTermType::StressForm, "Normal traction BC only supported for stress form of viscous term");
      zero_mean = false;
      this->AddVelNeumannBC(f,attr,pressure);
   }

   // Add a traction bc, this is essentially AddVelNeumannBC with zero_mean = false
   void AddTractionBC(mfem::VectorCoefficient *coeff, mfem::Array<int> &attr)
   {
      MFEM_ASSERT(viscous_term == ifs_solver::ViscousTermType::StressForm, "Traction BC only supported for stress form of viscous term");
      zero_mean = false;
      this->AddVelNeumannBC(coeff,attr);
   }
   void AddTractionBC(VecFuncT *f, mfem::Array<int> &attr)
   {
      MFEM_ASSERT(viscous_term == ifs_solver::ViscousTermType::StressForm, "Traction BC only supported for stress form of viscous term");
      zero_mean = false;
      this->AddVelNeumannBC(f,attr);
   }
   void AddTractionFreeBC(mfem::Array<int> &attr)
   {
      if (verbose && my_rank == 0)
      {
         mfem::out << "Adding Traction Free BC to attributes ";
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
         MFEM_ASSERT((vel_nat_attr[i] && attr[i]) == 0,
                     "Duplicate boundary definition detected.");
         if (attr[i] == 1)
         {
            vel_nat_attr[i] = 1;
         }
      }
      // do nothing condition
      zero_mean = false;
   }
   // Add outflow/homogeneous neumann boundary condition
   void AddOutflowBC(mfem::Array<int> &attr);

   //////////////// Pressure BCs ////////////////
   /// Add a Dirichlet boundary condition to the pressure field.
   void AddPresDirichletBC(mfem::Coefficient *coeff, mfem::Array<int> &attr);
   void AddPresDirichletBC(ScalarFuncT *f, mfem::Array<int> &attr);

   /// Set a nodal value of pressure to a constant
   /*
      Note this function expects coords to correspond to the 
      position of a node in the mesh!
   */
   void AddPressureNodeCondition(const double &p, mfem::Vector &coords, double tol = 1e-12);
   void AddPressureNodeCondition(ScalarFuncT &p_coeff, mfem::Vector &coords, double tol = 1e-12);

   // Enforce a zero mean pressure condition
   void AddZeroMeanPressureCondition()
   {
      if (verbose && my_rank == 0) mfem::out << "Adding zero mean pressure condition" << std::endl;
      zero_mean = true;
   }

   //////////////// Mixed BCs or Alias Names ////////////////
   void AddPressureBC(mfem::Coefficient *coeff, mfem::Array<int> &attr)
   {
      AddNormalTractionBC(coeff,attr,true);
      // do nothing in other components
   };
   
   //////////////// Body Forces ////////////////
   /// Add an body force to the RHS of the equation.
   /**
    * The VecFuncT @a f is evaluated at the current time t and extrapolated
    * together with the nonlinear parts of the Navier Stokes equation.
    */
   void AddBodyForce(mfem::VectorCoefficient *coeff, mfem::Array<int> &attr);
   void AddBodyForce(VecFuncT *f, mfem::Array<int> &attr);
   std::vector<AccelTerm_T> *GetBodyForces()
   {
      return &accel_terms;
   }

   ////////////////////////// Cardiovascular BCs and Related Functions //////////////////////////
   // backflow stabilization option
   void UseBackflowStabilization(mfem::Array<int> &attr, double beta_ = 0.2);

   // windkessel bc
   #ifdef MFEM_USE_MPI
   void AddWindkesselBC(mfem::Array<int> &attr, double R_, double Rd_, double C_, mfem::Coefficient *p_0_, std::function<double(double)> p_d_, double init_flow_rate_, mfem::ParFiniteElementSpace *gfes = nullptr);
   void AddWindkesselBC(mfem::Array<int> &attr, double R_, double Rd_, double C_, double p_0_, std::function<double(double)> p_d_, double init_flow_rate_, mfem::ParFiniteElementSpace *gfes = nullptr);
   void AddWindkesselBC(mfem::Array<int> &attr, double R_, double Rd_, double C_, double p_0_, double p_d_, double init_flow_rate_, mfem::ParFiniteElementSpace *gfes = nullptr);
   #else
   void AddWindkesselBC(mfem::Array<int> &attr, double R_, double Rd_, double C_, mfem::Coefficient *p_0_, std::function<double(double)> p_d_, double init_flow_rate_, mfem::FiniteElementSpace *gfes = nullptr);
   void AddWindkesselBC(mfem::Array<int> &attr, double R_, double Rd_, double C_, double p_0_, std::function<double(double)> p_d_, double init_flow_rate_, mfem::FiniteElementSpace *gfes = nullptr);
   void AddWindkesselBC(mfem::Array<int> &attr, double R_, double Rd_, double C_, double p_0_, double p_d_, double init_flow_rate_, mfem::FiniteElementSpace *gfes = nullptr);
   #endif

   // resistance bc
   #ifdef MFEM_USE_MPI
   void AddResistanceBC(mfem::Array<int> &attr, double R_, double p0_, mfem::ParFiniteElementSpace *gfes = nullptr);
   #else
   void AddResistanceBC(mfem::Array<int> &attr, double R_, double p0_, mfem::FiniteElementSpace *gfes = nullptr);
   #endif   

   #ifdef MFEM_USE_MPI
   // set mapped inflow bc
   enum mapped_inflow_type {Parabolic = 0, Womersley = 1}; // Womersley TO DO!
   MappedInflowCoefficient *mapped_inflow = nullptr;
   // Add mapped inflow bc boundary condition
   int inflow_type;
   // Note: this function requires the full mesh, so DO NOT delete mesh before setting this BC
   // time depend refers to radial inflow profile, aka the set inflow type
   void AddMappedInflowBC(int inflow_type_, Mesh *mesh, std::function<double (double)> flow_rate, bool time_depend, Array<int> &attr);
   // same function but using flow rate from file
   void AddMappedInflowBC_FromFile(int inflow_type_, Mesh *mesh, std::string filename_, bool time_depend, Array<int> &attr, int approx_type, bool periodic_data, int num_fourier_modes_ = 0);
   #endif

   ////////////////////////// Immersed Solid and Related Functions //////////////////////////
   #ifdef MFEM_USE_MPI
   void AddImmersedSolid(mfem::ParMesh *solid_pmesh, int orderd, mfem::Coefficient *init_rhos, MaterialType *mat_, int form_type_, bool subtrac_visc_stress = true, double scale_fac_ = 1.0);
   #else
   void AddImmersedSolid(mfem::Mesh *solid_pmesh, int orderd, mfem::Coefficient *init_rhos, MaterialType *mat_, int form_type_, bool subtrac_visc_stress = true, double scale_fac_ = 1.0);
   #endif
   void AddImmersedSolid(ImmersedSolid *im_solid_, double scale_fac_ = 1.0);

   //////////////////////
   /// SOLVER OPTIONS ///
   //////////////////////

   /// ADAPTIVE TIME STEPPING
   // set function, with option to set CFL num
   // est dt function
   // method for adaptive time stepping based on CFL_Crit
   double CFL_crit = 1.0;
   void SetCFL_Crit(double crit) {CFL_crit = crit;};
   double estimate_convective_cfl_dt();
   // TO DO^, code can handle changing dt, but need a way to choose dt   

   // set rho_infty for gen alpha
   void SetRhoInf(double rhoinf_ = 0.5)
   { 
      rho_inf = rhoinf_;
   }

   // flag for viscous term form
   int viscous_term = ifs_solver::ViscousTermType::StressForm;
   void SetViscousTermType (const int &type)
   {
      viscous_term = type;
   }

   // flag for dimensional form of equations
   int dimensional_flag = ifs_solver::dim_flag::NonDimensional;
   
   // Set nonlinear form type for CG
   void SetNonlinearFormType(const int& type)
   {
      nonlinearformtype = type;
   };
   int nonlinearformtype = ifs_solver::NonlinearFormType::Convective; // set default value

   // use a direct solver for the linear solve
   void UseDirectSolver()
   {
      direct_solver = true;
   };
   bool direct_solver = false;

   // Allocate a direct solver
   virtual void AllocDirectSolver();

   // Set preconditioner form type for CG
   void SetPreconditioner(const int& type, bool use_direct_solver = false)
   {
      prec_type = type;
      direct_solver_prec = use_direct_solver;
   };
   bool direct_solver_prec = false;
   int prec_type = ifs_solver::PreconditionerType::SIMPLE; // set default value

   // Set time stepping method
   void SetTimeSteppingMethod(const int& type)
   {
      time_stepping_method = type;
   };
   int time_stepping_method = ifs_solver::TimeSteppingMethod::GeneralizedAlpha; // set default value

   // Set stabilization type
   void SetStabilizationType(const int& type)
   {
      stab_type = type;
      switch (stab_type) {
         case ifs_solver::StabilizationType::NONE:
            break;

		   case ifs_solver::StabilizationType::SUPG:
            if (orderv != orderp - 1) 
            {
                  MFEM_ABORT("Error: using SUPG only requires mixed order elements");
            }
            this->UseSUPG();
            break;

		   case ifs_solver::StabilizationType::SUPG_PSPG:
            if (orderv != orderp) 
            {
                  MFEM_ABORT("Error: SUPG/PSPG requires equal order elements");
            }
            this->UseSUPG();
            this->UsePSPG();
            break;

		   case ifs_solver::StabilizationType::VMS:
         default:
            if (orderv != orderp) 
            {
                  MFEM_ABORT("Error: VMS requires equal order elements");
            }
            this->UseVMS();
            break;
      }
   };
   int stab_type = ifs_solver::StabilizationType::NONE;

   // options for stabilization
   void UseSUPG()
   {
      supg_flag = true;
   }

   void UsePSPG()
   {
      MFEM_ASSERT(equal_order,"must be equal order elements");
      pspg_flag = true;
   }

   void UseVMS()
   {
      MFEM_ASSERT(equal_order,"must be equal order elements");
      vms_flag = true;

      // also turn on PSPG and SUPG and mass stab?
      supg_flag = true;
      pspg_flag = true;
      lsic_flag = true;
   }

   void UseMassStabilization()
   {
      MFEM_ASSERT(supg_flag, "only supported with at least SUPG!");
      lsic_flag = true;
   }

   void UseTotalIntegrator()
   {
      total_integ = true;
      // use integrator that contains all terms in one!
   }

   // solver settings
   // set tolerances
   void SetLinearSolveAbsTol(double atol)
   {
      lin_atol = atol;
   };
   void SetLinearSolveRelTol(double rtol)
   {
      lin_rtol = rtol;
   };
   void SetLinearSolveMaxIt(int maxIt)
   {
      lin_maxIt = maxIt;
   };
   void SetLinearSolveMinIt(int minIt)
   {
      MFEM_ABORT("TO DO");
      lin_minIt = minIt;
   };
   void SetNewtonSolveAbsTol(double atol)
   {
      atol_n = atol;
   };
   void SetNewtonSolveRelTol(double rtol)
   {
      rtol_n = rtol;
   };
   void SetNewtonSolveMaxIt(int maxIt)
   {
      maxIt_n = maxIt;
   };
   void SetNewtonSolveMinIt(int minIt)
   {
      minIt_n = minIt;
   };  
   void SetNewtonSolveNoErrorNonConvergence()
   {
      ns_no_error_nonconvergence = true;
   };

   void NoPreconditioner()
   {
      no_prec = true;
   };
   
   void UseSIMPLEPreconditioner()
   {
      simple_prec_flag = true;
   };

   void UseSCRPreconditioner()
   {
      scr_prec_flag = true;
   };

   void Use_aSCRPreconditioner()
   {
      ascr_prec_flag = true;
   };

   void PrintNewton(bool val)
   {
      print_newton = val;
   };

   // ainv and sinv on preconditioners
   void SetPrecAINVAbsTol(double atol)
   {
      ainv_prec_atol = atol;
   };
   void SetPrecAINVRelTol(double rtol)
   {
      ainv_prec_rtol = rtol;
   };

   void SetPrecAINV_Inner_AbsTol(double atol)
   {
      ainv_inner_atol = atol;
   };
   void SetPrecAINV_Inner_RelTol(double rtol)
   {
      ainv_inner_rtol = rtol;
   };

   void SetPrecSINVAbsTol(double atol)
   {
      sinv_prec_atol = atol;
   };
   void SetPrecSINVRelTol(double rtol)
   {
      sinv_prec_rtol = rtol;
   };


   ///////////////////////
   /// POST-PROCESSING ///
   ///////////////////////

   /// Compute \f$\nabla \times \nabla \times u\f$ for \f$u \in (H^1)^2\f$.
   void ComputeCurl2D(mfem::GridFunction &u,
                      mfem::GridFunction &cu,
                      bool assume_scalar = false);

   /// Compute \f$\nabla \times \nabla \times u\f$ for \f$u \in (H^1)^3\f$.
   void ComputeCurl3D(mfem::GridFunction &u, mfem::GridFunction &cu);

   // compute wall-shear stress vector
   /// Return a pointer to the current wall shear stress vector ParGridFunction.
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction *GetCurrentWSS() { return &wss_gf; }
   #else
   mfem::GridFunction *GetCurrentWSS() { return &wss_gf; }
   #endif
   // Computes WSS using either simple averaging (Standard) or Lumped L2 Projection (Row-sum mass lumping)
   void ComputeWSS(mfem::GridFunction &u, mfem::GridFunction &p, int avg_type = ifs_solver::WSSAvgType::Standard);

   /// Compute CFL
   double ComputeCFL(const mfem::GridFunction &u, double dt);

   // Compute integral of div v over domain
   double ComputeIntegralDivV(const mfem::GridFunction &v_gf);

   double ComputeMassFlowThruSurf(mfem::GridFunction &u, mfem::Array<int> &attr_surf, bool inflow_flag, double fluid_density, bool volume_flow_rate = false);

   double ComputeAvgPressureOnSurf(mfem::GridFunction &p_gf, mfem::Array<int> &attr_surf, double &area);

   double ComputeAvgVeloOnSurf(mfem::GridFunction &u, mfem::Array<int> &attr_surf, mfem::Vector &avg_velo);

   //////////////////
   //// RESTART /////
   //////////////////
   void Restarting(std::string restart_folder_, int step, double &dt, double &time, bool ics = false);
   void writeRestartFiles(std::string restart_folder_, int restart_step_size_);
   void writeRestarts(int step, double dt, double time);
   void readRestarts(int &step, double &dt, double &time, bool read_bcs = true);

   //////////////////
   /// DESTRUCTOR ///
   //////////////////

   virtual ~IncompressibleFluidSolver();

protected:
   /// Print information about the Navier version.
   void PrintInfo();

   ////////////////
   /// SETTINGS ///
   ////////////////

   // choice of ordering solution unknowns
   // int orderingtype = Ordering::byVDIM;
   int orderingtype = mfem::Ordering::byNODES;

   /// Enable/disable debug output.
   bool debug = false;

   /// Enable/disable verbose output.
   bool verbose = true;

   // zero mean
   // default is false so a traction free condition is applied when no bc is specified
   bool zero_mean = false; 
   double ref_p = 0.0;
   int p_local_node_num = -1;
   bool pressure_node_condition = false;
   ScalarFuncT *ref_p_coeff = nullptr; // for time dependent pressure node value // not owned

   // outflow_bc_flag
   bool outflow_bc_flag = false;
   mfem::Array<int> outflow_attr;

   // Stabilization options
   // enable/disable equal order interpolation
   bool equal_order = false;
   bool supg_flag = false;
   bool pspg_flag = false;
   bool vms_flag = false;
   bool lsic_flag = false;

   /// The mesh.
   #ifdef MFEM_USE_MPI
   mfem::ParMesh *pmesh = nullptr; // not owned
   #else
   mfem::Mesh *pmesh = nullptr; // not owned
   #endif
   int my_rank = 0;

   /// The order of the velocity and pressure space.
   int orderv;
   int orderp;

   // vector dimension for velocity field
   int vdim;

   /// Kinematic viscosity (dimensionless).
   double kin_vis;
   double fluid_density;
   mfem::ConstantCoefficient rhof; // only support constant fluid density for now   
   mfem::ConstantCoefficient mu_f;

   /// Velocity \f$H^1\f$ finite element collection.
   mfem::FiniteElementCollection *vfec = nullptr;

   /// Pressure \f$H^1\f$ finite element collection.
   mfem::FiniteElementCollection *pfec = nullptr;

   /// Velocity \f$(H^1)^d\f$ finite element space.
   #ifdef MFEM_USE_MPI
   mfem::ParFiniteElementSpace *vfes = nullptr;
   #else
   mfem::FiniteElementSpace *vfes = nullptr;
   #endif

   /// Pressure \f$H^1\f$ finite element space.
   #ifdef MFEM_USE_MPI
   mfem::ParFiniteElementSpace *pfes = nullptr;
   #else
   mfem::FiniteElementSpace *pfes = nullptr;
   #endif

   /// Linear form to compute the mass matrix in various subroutines.
   #ifdef MFEM_USE_MPI
   mfem::ParLinearForm *mass_lf = nullptr;
   #else
   mfem::LinearForm *mass_lf = nullptr;
   #endif
   mfem::ConstantCoefficient onecoeff;
   double volume = 0.0;

   // All essential attributes
   mfem::Array<int> vel_ess_attr;
   mfem::Array<int> pres_ess_attr;

   // All natural attributes
   mfem::Array<int> vel_nat_attr;

   // Periodic attributes
   mfem::Array<int> periodic_attr;

   // Neumann attributes
   mfem::Array<int> neumann_attr;

   // Bookkeeping for velocity dirichlet bcs.
   std::vector<VelDirichletBC_T> vel_dbcs;

   // Bookkeeping for velocity neumann bcs.
   std::vector<VelNeumannBC_T> vel_nbcs;

   // Bookkeeping for pressure dirichlet bcs.
   std::vector<PresDirichletBC_T> pres_dbcs;

   // Bookkeeping for acceleration (forcing) terms.
   std::vector<AccelTerm_T> accel_terms;

   // essential dof arrays
   mfem::Array<int> vel_ess_tdofs;
   mfem::Array<int> pres_ess_tdofs;

   // backflow stabilization option
   bool backflow_stab = false;
   // bookkeeping 
   std::vector<BackflowStabilization_T> backflow_stabilization_bcs;

   ///// Cardiovascular BCs /////
   // resistance bc option
   bool resistance_bc = false;
   // Bookkeeping for resistance bcs
   std::vector<ResistanceBC_T> resistance_bcs;

   // windkessel
   bool windkessel_bc = false;
   // Bookkeeping for windkessel bcs
   std::vector<WindkesselBC_T> windkessel_bcs;

   CVBC_Op *cvbc_op = nullptr;
   mfem::Vector A3d_diag, Abc_diag, Atotal_diag;
   CVBC_AinvPC *cvbc_PC = nullptr;

   // ALE
   bool ale = false;

   ///// Immersed Solid /////
   bool IFEM = false;
   immersed_solid::ImmersedSolid* im_solid = nullptr;
   mfem::Vector f_fsi;
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction s_gf;
   #else
   mfem::GridFunction s_gf;
   #endif
   double scale_fac = 1.0;

   // Time stepper data
   int cur_step = 0;
   double curr_dt = 0.0;

   // Generalized alpha data
   double alpha_m = 0.0;
   double alpha_f = 0.0;
   double gamma = 0.0;
   double beta = 0.0;
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction vdot_gf;
   #else
   mfem::GridFunction vdot_gf;
   #endif
   double rho_inf = 0.5; // default
   
   // Spatial Operators
   tau_m_coeff tau_m;
   nu_LSIC nu_coeff;
   bool nxx_flag = false;
   #ifdef MFEM_USE_MPI
   mfem::Array<mfem::ParFiniteElementSpace *> fes_stab;
   #else
   mfem::Array<mfem::FiniteElementSpace *> fes_stab;
   #endif
   mfem::Array<mfem::Array<int> *> total_ess_tdofs;
   #ifdef MFEM_USE_MPI
   mfem::ParLinearForm *f_form = nullptr;
   mfem::ParLinearForm *L_L_form = nullptr;
   #else
   mfem::LinearForm *f_form = nullptr;
   mfem::LinearForm *L_L_form = nullptr;
   #endif
   ResidualOpFluid *res = nullptr;

   bool total_integ = true;
   #ifdef MFEM_USE_MPI
   mfem::ParBlockNonlinearForm *FluidNF = nullptr;
   #else
   mfem::BlockNonlinearForm *FluidNF = nullptr;
   #endif

   // Solvers and Preconditioners
   FluidNewtonSolver *newton_solver = nullptr; // newton solver
   double atol_n = 1e-10;
   double rtol_n = 1e-6;
   int maxIt_n = 10;
   int minIt_n = 0;
   int pl_n = -1;
   bool print_newton = true;
   bool ns_no_error_nonconvergence = false;

   // nonlinear residual final values
   double norm = 0.0;
   double norm_mass = 0.0;
   double norm_mom = 0.0;
   double rel_norm = 1.0;
   double rel_norm_mass = 1.0;
   double rel_norm_mom = 1.0;

   mfem::Solver *solver = nullptr; // linear solver
   double lin_atol = 1e-10;
   double lin_rtol = 1e-6;
   int lin_maxIt = 200;
   int lin_pl_solve = -1;
   int lin_minIt = 0;
   IterativeSolverMonitor *solver_monitor = nullptr;

   BlockOrthoSolver *solver_wrap = nullptr; // wrapper to apply mean zero condition

   mfem::Solver *jacobian_prec = nullptr; // block preconditioner
   bool simple_prec_flag = false;
   bool scr_prec_flag = false;
   bool ascr_prec_flag = false;
   bool no_prec = false;
   bool monolithic_prec = false;
   #ifdef MFEM_USE_MPI
   bool use_jacobi_prec = false; // option to use jacobi instead of AMG
   #else
   bool use_jacobi_prec = true;
   #endif
   // tolerances on inner solves for preconditioner
   double ainv_prec_atol = 1e-10;
   double ainv_prec_rtol = 1e-6;
   double sinv_prec_atol = 1e-10;
   double sinv_prec_rtol = 1e-6;

   double ainv_inner_atol = 1e-10; // only for scr!
   double ainv_inner_rtol = 1e-6;   

   // Solution Vectors and GridFunctions
   mfem::Vector ft, fb;
   
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction vn_gf, pn_gf, vdotn_gf, pdotn_gf;
   mfem::ParGridFunction vn_next_gf, pn_next_gf;
   #else
   mfem::GridFunction vn_gf, pn_gf, vdotn_gf, pdotn_gf;
   mfem::GridFunction vn_next_gf, pn_next_gf;
   #endif
   mfem::Vector temp_vec;

   // block solution vectors
   mfem::BlockVector yn, ydotn, ynp1, ydotnp1, yg, ydotg;
   mfem::Array<int> block_trueOffsets;
   mfem::Device device;
   mfem::MemoryType mt;

   // Post-Process
   // timing
   void RegisterTimers();
   TimerRegistry timer_reg;
   mfem::StopWatch sw_solve, sw_setup, sw_tan, sw_res, sw_linsolve, sw_initiators, sw_correctors;
   mfem::StopWatch sw_fluidtan, sw_fluidres, sw_solidtan, sw_solidres, sw_fsires, sw_fsitan, sw_fsicmat;
   mfem::StopWatch sw_predictors, sw_essbcs, sw_weakbcs, sw_restarts, sw_updates;
   bool extended_timing = true;
   bool timers_registered = false;

   // solver outputs
   // Iteration counts.
   int iter_solve = 0;

   // Residuals.
   double res_solve = 0.0;
   double rel_res_solve = 0.0;

   // GFs
   #ifdef MFEM_USE_MPI
   mfem::ParGridFunction wss_gf;
   #else
   mfem::GridFunction wss_gf;
   #endif

   // restarting data
   std::string restart_folder = "RestartFiles";
   int restart_save_step = 1;
   int restart_step = 1;
   bool write_restarts = false;
   bool restarting = false;
};

} // namespace ins_solver_serial