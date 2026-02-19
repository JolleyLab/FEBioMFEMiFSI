#include "IncompressibleFluidSolver.hpp"
#include <general/forall.hpp>
#include <fstream>
#include <iomanip>
#include "../Utilities/linear_element_ext.hpp"
#include "../Utilities/CPardisoSolverExt.hpp"
#include <memory>
#include <filesystem>
#include <cmath>

using namespace mfem;
using namespace ifs_solver;
using namespace std;
using namespace immersed_solid;

#ifdef MFEM_USE_MPI
IncompressibleFluidSolver::IncompressibleFluidSolver(ParMesh *mesh, int orderv, int orderp, double kin_vis, int ordering_type_, bool immersed_flag)
   : pmesh(mesh), orderv(orderv), orderp(orderp), kin_vis(kin_vis), orderingtype(ordering_type_)
#else
IncompressibleFluidSolver::IncompressibleFluidSolver(Mesh *mesh, int orderv, int orderp, double kin_vis, int ordering_type_, bool immersed_flag)
   : pmesh(mesh), orderv(orderv), orderp(orderp), kin_vis(kin_vis), orderingtype(ordering_type_)
#endif
{
   // select backend device   
   // TODO: this should probably be passed into the constructor and set outside the class
   // also any vectors should be allocated with SetSize(size, memtype)!
   const char *device_config = "cpu"; // change as needed!
   // #ifdef MFEM_USE_OPENMP
   // const char *device_config = "omp";
   // #else
   // const char *device_config = "cpu";
   // #endif
   device.Configure(device_config);
   mt = device.GetMemoryType();

   // check type of elements in mesh
   Array<Geometry::Type> geom_types;
   mesh->GetGeometries(pmesh->Dimension(), geom_types);
   if (pmesh->Dimension() == 3)
   {
      for (int i = 0; i < geom_types.Size(); ++i)
      {
         if (geom_types[i] == Geometry::PYRAMID || geom_types[i] == Geometry::PRISM)
         {
            nxx_flag = false;
         }
      }
   }
   else
   {
      nxx_flag = true;
   }

   MFEM_ASSERT(orderp >= 1, "Order of pressure approx must be >= 1");
   equal_order = orderv == orderp ? true : false;
   vfec = new H1_FECollection(orderv, pmesh->Dimension());
   pfec = new H1_FECollection(orderp, pmesh->Dimension()); 
   #ifndef MFEM_USE_MPI
   vfes = new FiniteElementSpace(pmesh, vfec, pmesh->Dimension(),orderingtype);
   pfes = new FiniteElementSpace(pmesh, pfec,1,orderingtype);
   #else
   vfes = new ParFiniteElementSpace(pmesh, vfec, pmesh->Dimension(),orderingtype);
   pfes = new ParFiniteElementSpace(pmesh, pfec,1,orderingtype);
   my_rank = pmesh->GetMyRank();
   #endif

   vdim = pmesh->Dimension();

   // for isoparametric mapping
   pmesh->SetNodalFESpace(vfes); 

   // Check if fully periodic mesh
   if (!(pmesh->bdr_attributes.Size() == 0))
   {
      // dirichlet bcs
      vel_ess_attr.SetSize(pmesh->bdr_attributes.Max() * vdim);
      vel_ess_attr = 0;
      
      pres_ess_attr.SetSize(pmesh->bdr_attributes.Max());
      pres_ess_attr = 0;

      // neumann bcs
      vel_nat_attr.SetSize(pmesh->bdr_attributes.Max());
      vel_nat_attr = 0;

      neumann_attr.SetSize(pmesh->bdr_attributes.Max());
      neumann_attr = 0;
      
      // periodic bcs
      periodic_attr.SetSize(pmesh->bdr_attributes.Max());
      periodic_attr = 0;

      outflow_attr.SetSize(pmesh->bdr_attributes.Max());
      outflow_attr = 0;
   }
   else
   {
      periodic_attr.SetSize(pmesh->bdr_attributes.Max());
      periodic_attr = 1;
   }

   int vfes_truevsize = vfes->GetTrueVSize();
   int pfes_truevsize = pfes->GetTrueVSize();

   ft.SetSize(vfes_truevsize);
   ft = 0.0;
   fb.SetSize(vfes_truevsize);
   fb = 0.0;

   vn_gf.SetSpace(vfes);
   vn_gf = 0.0;
   vdotn_gf.SetSpace(vfes);
   vdotn_gf = 0.0;
   pn_gf.SetSpace(pfes);
   pn_gf = 0.0;
   pdotn_gf.SetSpace(pfes);
   pdotn_gf = 0.0;
   vn_next_gf.SetSpace(vfes);
   vn_next_gf = 0.0;
   pn_next_gf.SetSpace(pfes);
   pn_next_gf = 0.0;

   // block vectors
   block_trueOffsets.SetSize(3); // number of variables + 1
   block_trueOffsets[0] = 0;
   block_trueOffsets[1] = vfes->GetTrueVSize();
   block_trueOffsets[2] = pfes->GetTrueVSize();
   block_trueOffsets.PartialSum();

   // block vectors containing yn = [un pn] and ydotn = [udotn pdotn]
   yn.Update(block_trueOffsets,mt);
   yn = 0.0;
   ynp1.Update(block_trueOffsets,mt);
   ynp1 = 0.0;
   ydotn.Update(block_trueOffsets,mt);
   ydotn = 0.0;
   ydotnp1.Update(block_trueOffsets,mt);
   ydotnp1 = 0.0;
   yg.Update(block_trueOffsets,mt);
   yg = 0.0;
   ydotg.Update(block_trueOffsets,mt);
   ydotg = 0.0;

   cur_step = 0;

   PrintInfo();

   fluid_density = 1.0;
   rhof.constant = fluid_density;
}

void IncompressibleFluidSolver::MeanZero(GridFunction &v)
{  
   // Make sure not to recompute the inner product linear form every
   // application. Only update if mesh moves!
   if (mass_lf == nullptr)
   {
      onecoeff.constant = 1.0;
      #ifndef MFEM_USE_MPI
      mass_lf = new LinearForm(v.FESpace());
      #else
      mass_lf = new ParLinearForm(dynamic_cast<ParGridFunction*>(&v)->ParFESpace());
      #endif
      mass_lf->AddDomainIntegrator(new DomainLFIntegrator(onecoeff));
      mass_lf->Assemble();

      #ifndef MFEM_USE_MPI
      GridFunction one_gf(v.FESpace());
      #else
      ParGridFunction one_gf(dynamic_cast<ParGridFunction*>(&v)->ParFESpace());
      #endif
      one_gf.ProjectCoefficient(onecoeff);

      volume = mass_lf->operator()(one_gf);
   }

   #ifdef MFEM_USE_MPI
   double integ = mass_lf->operator()(*dynamic_cast<ParGridFunction*>(&v));
   #else
   double integ = mass_lf->operator()(v);
   #endif

   v -= integ / volume;
}

void IncompressibleFluidSolver::Orthogonalize(Vector &v)
{
   double global_sum = v.Sum();
   int global_size = v.Size();

   #ifdef MFEM_USE_MPI
   MPI_Allreduce(MPI_IN_PLACE, &global_sum, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &global_size, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   #endif

   v -= global_sum / static_cast<double>(global_size);
}

void IncompressibleFluidSolver::AllocDirectSolver()
{
   #ifndef MFEM_USE_MPI
      #ifdef MFEM_USE_SUITESPARSE
         solver = new UMFPackSolver();
         dynamic_cast<UMFPackSolver*>(solver)->Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      #else
         MFEM_ABORT("Error: Serial direct solvers require mfem with SuiteSparse!");
      #endif
   #else
      #ifdef MFEM_USE_MKL_CPARDISO
         if (IFEM) 
         {
            solver = new CPardisoSolverExt(vfes->GetComm());
            dynamic_cast<CPardisoSolverExt*>(solver)->SetMatrixType(CPardisoSolverExt::MatType::REAL_NONSYMMETRIC);
            dynamic_cast<CPardisoSolverExt*>(solver)->ReinitFlag();
         }
         else
         {
            solver = new CPardisoSolver(vfes->GetComm());
            dynamic_cast<CPardisoSolver*>(solver)->SetMatrixType(CPardisoSolver::MatType::REAL_NONSYMMETRIC);
         }   
      #else
         MFEM_ABORT("Error: Parallel direct solvers require mfem with MKL PARDISO!");
      #endif
   #endif
}

void IncompressibleFluidSolver::Setup(double dt)
{
   sw_setup.Start();
   Array<int> empty; // for forming some of the matrices

   // remove duplicated dofs from appends!
   vel_ess_tdofs.Unique();
   pres_ess_tdofs.Unique();

   vel_ess_tdofs.Sort();
   pres_ess_tdofs.Sort();

   /////////////////////////////////////
   ////// Setup Spatial Operators //////
   /////////////////////////////////////

   // set entire neumann boundary
   for (int i = 0; i < neumann_attr.Size(); ++i)
   {
      // check for specified neumann bcs
      if (pres_ess_attr[i] == 1 || vel_nat_attr[i] == 1 || outflow_attr[i] == 1)
      {
         neumann_attr[i] = 1;
      }

      // check resistance bcs
      if (resistance_bc)
      {
         for (auto &rbc : resistance_bcs)
         {
            if (rbc.attr[i] == 1)
            {
               neumann_attr[i] = 1;
            }
         }  
      }

      // check windkessel bcs
      if (windkessel_bc)
      {
         for (auto &wbc : windkessel_bcs)
         {
            if (wbc.attr[i] == 1)
            {
               neumann_attr[i] = 1;
            }
         }  
      }
   }

   // set gen alpha parameters
   if (time_stepping_method == ifs_solver::TimeSteppingMethod::GeneralizedAlpha)
   {  
      alpha_m = 0.5*( (3.0 - rho_inf) / (1.0 + rho_inf) );
      alpha_f = 1.0 / (1.0 + rho_inf);
      gamma = 0.5 + alpha_m - alpha_f;
      beta = 0.25 * (1.0 - alpha_f + alpha_m) * (1.0 - alpha_f + alpha_m); // for immersed fsi, so only one method is used!
   }
   else if (time_stepping_method == ifs_solver::TimeSteppingMethod::BackwardEuler)
   {
      // BE from Gen-Alpha
      alpha_m = 1.0;
      alpha_f = 1.0;
      gamma = 1.0;
      beta = 9.0/16.0;
   }
   else
   {
      MFEM_ABORT("Error: Unknown time stepping method");
   }   

   if (verbose && my_rank == 0)
   {
      mfem::out << "Setup fluid" << std::endl;
      if (direct_solver) mfem::out << "SOLVER: Direct Solver" << std::endl;
      else
      {
         mfem::out << "SOLVER: FGMRES" << std::endl;
         if (prec_type == ifs_solver::PreconditionerType::SIMPLE) mfem::out << "Preconditioner = SIMPLE" << std::endl;
         if (prec_type == ifs_solver::PreconditionerType::SCR) mfem::out << "Preconditioner = SCR" << std::endl;
         if (prec_type == ifs_solver::PreconditionerType::aSCR) mfem::out << "Preconditioner = a-SCR" << std::endl;
         if (prec_type == ifs_solver::PreconditionerType::ILU) mfem::out << "Preconditioner = ILU" << std::endl;
         if (no_prec) mfem::out << "Preconditioner = NONE" << std::endl;
      }
      if (nonlinearformtype == ifs_solver::NonlinearFormType::Convective) mfem::out << "Nonlinear term TYPE = Convective" << std::endl;
      if (nonlinearformtype == ifs_solver::NonlinearFormType::EMAC) mfem::out << "Nonlinear term TYPE = EMAC" << std::endl;
      if (viscous_term == ifs_solver::ViscousTermType::LaplaceForm) mfem::out << "Viscous term TYPE = LaplaceForm" << std::endl;
      if (viscous_term == ifs_solver::ViscousTermType::StressForm) mfem::out << "Viscous term TYPE = StressForm" << std::endl;
      else mfem::out << "Enforce avg p = 0 : FALSE" << std::endl;
      if (supg_flag || pspg_flag || lsic_flag || vms_flag) mfem::out << "Stabilizations:" << std::endl;
      if(supg_flag) mfem::out << "SUPG" << std::endl;
      if(pspg_flag) mfem::out << "PSPG" << std::endl;
      if(vms_flag) mfem::out << "VMS" << std::endl;
      if(lsic_flag) mfem::out << "LSIC" << std::endl;
      
      mfem::out << "TIME STEPPING = Generalized Alpha" << std::endl;
      mfem::out << "rho_inf = " << rho_inf << std::endl;
      mfem::out << "gamma = " << gamma << std::endl;
      mfem::out << "beta = " << beta << std::endl;
      mfem::out << "alpha_f = " << alpha_f << std::endl;
      mfem::out << "alpha_m = " << alpha_m << std::endl;
      mfem::out << "---------------------------------------" << std::endl;
   }

   if (total_integ)
   {
      mu_f.constant = dimensional_flag == ifs_solver::dim_flag::NonDimensional ? kin_vis : fluid_density*kin_vis;


      if (equal_order) MFEM_ASSERT(pspg_flag, "Error: Must use PSPG with equal order elements!");

      // tau_m = tau_PSPG = tau_SUPG
      tau_m.kinvis = kin_vis;
      tau_m.dt = &curr_dt;
      tau_m.Set_CI(orderv);
      if (vms_flag) 
      {
         tau_m.UseAdditionalMomentumRestoringTerm();
         // tau_m.SetTol(atol_n);
      }

      if (IFEM && s_gf.FESpace() != nullptr)
      {
         tau_m.s_gf = &s_gf;
         tau_m.scale_fac = scale_fac;
      }

      // Need to update this every newton iteration!
      vdot_gf.SetSpace(vfes);
      vdot_gf = 0.0;

      if (lsic_flag)
      {
         nu_coeff.tau_coeff = &tau_m;
      }

      fes_stab.SetSize(2);
      fes_stab[0] = vfes;
      fes_stab[1] = pfes;

      auto IncFluidInteg = new IncompressibleFluidIntegrator();
      IncFluidInteg->SetFluidProperties(rhof, &mu_f);
      IncFluidInteg->SetFluidData(&vdot_gf, &accel_terms, viscous_term);
      IncFluidInteg->SetTimeData(alpha_m, alpha_f, gamma, &curr_dt);
      IncFluidInteg->SetStabilizationTerms(supg_flag, pspg_flag, vms_flag, lsic_flag, tau_m, nu_coeff);
      IncFluidInteg->NxxFlag(nxx_flag);
      #ifndef MFEM_USE_MPI
      FluidNF = new BlockNonlinearForm(fes_stab);
      #else
      FluidNF = new ParBlockNonlinearForm(fes_stab);
      #endif
      FluidNF->AddDomainIntegrator(IncFluidInteg); // FluidNF assumes ownership of IncFluidInteg
   }
   else
   {
      MFEM_ABORT("TO DO");
   }

   /// Surface tractions ///
   #ifndef MFEM_USE_MPI
   L_L_form = new LinearForm(vfes);
   #else
   L_L_form = new ParLinearForm(vfes);
   #endif
   // neumann velocity bcs/tractions
   for (auto &vel_nbc : vel_nbcs)
   {
      if(vel_nbc.coeff)
      {
         L_L_form->AddBoundaryIntegrator(new TractionIntegrator(*vel_nbc.coeff),vel_nbc.attr);
      }
      else if(vel_nbc.mcoeff)
      {
         L_L_form->AddBoundaryIntegrator(new TractionIntegrator(*vel_nbc.mcoeff),vel_nbc.attr);
      }
      else if (vel_nbc.scalar_coeff)
      {
         L_L_form->AddBoundaryIntegrator(new TractionIntegrator(*vel_nbc.scalar_coeff, vel_nbc.pressure), vel_nbc.attr);
      }
      else
      {
         MFEM_ABORT("Error in defining neumann bcs for velocity");
      }
   }

   // backflow stabilization
   if (backflow_stab)
   {
      MFEM_ASSERT(total_integ, "Error: Backflow stabilization not implemented without total_integ = true");
      for (auto &bf_bc : backflow_stabilization_bcs)
      {
         FluidNF->AddBoundaryIntegrator(new BlockBackflowStabIntegrator(&curr_dt,alpha_f,gamma,rhof,bf_bc.beta,true),bf_bc.attr);
      }   
   }  

   // cardiovascular bcs
   if (windkessel_bc || resistance_bc)
   {
      cvbc_op = new CVBC_Op(vfes->GetTrueVSize(),vfes);
      
      // resistance bcs
      if(resistance_bc) cvbc_op->SetResistanceBC(&resistance_bcs);

      // windkessel bcs
      if(windkessel_bc)
      {
         // initialize windkessel bcs
         for (auto &wbc : windkessel_bcs)
         {            
            wbc.q_n = wbc.init_flow_rate;  
            wbc.t_n = 0.0;
            wbc.dt = dt;  
            wbc.t_np1 = dt; 
            wbc.integ_hist_terms = 0.0;
         }  

         cvbc_op->SetWindkesselBC(&windkessel_bcs);
      }

      // vectors to store relevant diagonals for A block preconditioner with CVBCs
      A3d_diag.SetSize(vfes->GetTrueVSize());
      A3d_diag = 0.0;

      Abc_diag.SetSize(vfes->GetTrueVSize());
      Abc_diag = 0.0;

      Atotal_diag.SetSize(vfes->GetTrueVSize());
      Atotal_diag = 0.0;

      // Prec for A block
      cvbc_PC = new CVBC_AinvPC(vfes->GetTrueVSize(),&A3d_diag,&Abc_diag,&Atotal_diag,vfes);
      if (resistance_bc) cvbc_PC->SetResistanceBC(&resistance_bcs);
      if (windkessel_bc) cvbc_PC->SetWindkesselBC(&windkessel_bcs);    
   } 

   /// Setup immersed solids ///
   if (IFEM)
   {
      im_solid->Setup(dt,alpha_m,alpha_f,gamma,beta); 
      im_solid->SetTimers(&sw_solidres,&sw_solidtan,&sw_fsires,&sw_fsitan,&sw_fsicmat);       
   }

   /////////////////////////////////////////////////////////
   ////// Setup Block Preconditioner and Newton Solve //////
   /////////////////////////////////////////////////////////

   // setup linear solver
   if (direct_solver)
   {
      AllocDirectSolver();
      no_prec = true;
   }
   else
   {
      // #ifndef MFEM_USE_MPI
      // solver = new ExtendedFGMRESSolver();
      // #else
      // solver = new ExtendedFGMRESSolver(vfes->GetComm());
      // #endif
      // static_cast<ExtendedFGMRESSolver *>(solver)->SetMinIter(lin_minIt);
      // static_cast<ExtendedFGMRESSolver *>(solver)->SetAbsTol(lin_atol);
      // static_cast<ExtendedFGMRESSolver *>(solver)->SetRelTol(lin_rtol);
      // static_cast<ExtendedFGMRESSolver *>(solver)->SetMaxIter(lin_maxIt);
      // static_cast<ExtendedFGMRESSolver *>(solver)->SetPrintLevel(lin_pl_solve);
      // // static_cast<ExtendedFGMRESSolver *>(solver)->SetKDim(lin_maxIt);

      #ifndef MFEM_USE_MPI
      solver = new FGMRESSolver();
      #else
      solver = new FGMRESSolver(vfes->GetComm());
      #endif
      static_cast<FGMRESSolver *>(solver)->SetAbsTol(lin_atol);
      static_cast<FGMRESSolver *>(solver)->SetRelTol(lin_rtol);
      static_cast<FGMRESSolver *>(solver)->SetMaxIter(lin_maxIt);
      static_cast<FGMRESSolver *>(solver)->SetPrintLevel(lin_pl_solve);
      // static_cast<FGMRESSolver *>(solver)->SetKDim(lin_maxIt);

      // setup linear solver monitor
      solver_monitor = new MinItMonitor(lin_minIt);
      static_cast<IterativeSolver *>(solver)->SetMonitor(*solver_monitor);
      auto iterative_solver = static_cast<IterativeSolver *>(solver);
      solver_monitor->SetIterativeSolver(*iterative_solver);

      // setup jacobian preconditioner
      Array<FiniteElementSpace *> fes;
      fes.SetSize(2);
      fes[0] = vfes;
      fes[1] = pfes;

      // monolithic preconditioners
      if (prec_type == ifs_solver::PreconditionerType::Jacobi)
      {
         #ifdef MFEM_USE_MPI
         jacobian_prec = new HypreSmoother();
         dynamic_cast<HypreSmoother*>(jacobian_prec)->SetType(HypreSmoother::Jacobi, 1);

         MFEM_ASSERT(windkessel_bc == false && resistance_bc == false,
                     "Error: Jacobi preconditioner does not support CVBCs currently!");

         if (!equal_order)
         {
            MFEM_WARNING("Warning: A monolithic Jacobi preconditioner may not work well with mixed order method (Pk-Pk-1)!");
         }

         monolithic_prec = true;
         #else
         jacobian_prec = new DSmoother();
         #endif
      }
      else if (prec_type == ifs_solver::PreconditionerType::ILU)
      {
         #ifdef MFEM_USE_MPI
         jacobian_prec = new HypreILU();
         // jacobian_prec = new HypreEuclid(vfes->GetComm()); // Same as ILU but scalable

         MFEM_ASSERT(windkessel_bc == false && resistance_bc == false,
                     "Error: ILU preconditioner does not support CVBCs currently!");

         if (!equal_order)
         {
            MFEM_WARNING("Warning: A monolithic ILU preconditioner may not work well with mixed order method (Pk-Pk-1)!");
         }

         monolithic_prec = true;
         #else
         MFEM_ABORT("Error: ILU preconditioner requires mfem with MPI!");
         #endif
      }
      else
      {
         // block preconditioners
         if (prec_type == ifs_solver::PreconditionerType::SIMPLE)
         {
            jacobian_prec = new SIMPLE_Prec(block_trueOffsets,fes,equal_order);
         }
         else if (prec_type == ifs_solver::PreconditionerType::SCR)
         {
            #ifdef MFEM_USE_MPI
            jacobian_prec = new SCR_Prec(block_trueOffsets,fes,equal_order);
            dynamic_cast<SCR_Prec*>(jacobian_prec)->SetA_Inner_invTol(ainv_inner_atol,ainv_inner_rtol);
            #else
            MFEM_ABORT("TO DO");
            #endif         
         }
         else if (prec_type == ifs_solver::PreconditionerType::aSCR)
         {
            #ifdef MFEM_USE_MPI
            jacobian_prec = new aSCR_Prec(block_trueOffsets,fes,equal_order);
            #else
            MFEM_ABORT("TO DO");
            #endif 
         }
         else
         {
            MFEM_ABORT("Error: Preconditioner type not supported!");
         }

         // block preconditioner options
         if (windkessel_bc || resistance_bc) dynamic_cast<BlockPreconditioner*>(jacobian_prec)->SetCV_Ainv_PC(cvbc_PC);
         if (use_jacobi_prec) dynamic_cast<BlockPreconditioner*>(jacobian_prec)->UseJacobiPreconditioning();
         if (direct_solver_prec) dynamic_cast<BlockPreconditioner*>(jacobian_prec)->UseDirectSolvers();

         dynamic_cast<BlockPreconditioner*>(jacobian_prec)->SetAinvTol(ainv_prec_atol,ainv_prec_rtol);
         dynamic_cast<BlockPreconditioner*>(jacobian_prec)->SetSinvTol(sinv_prec_atol,sinv_prec_rtol);
      }      
   }
   
   // setup residual operator
   res = new ResidualOpFluid(block_trueOffsets,vel_ess_tdofs,pres_ess_tdofs,equal_order);
   res->SetTotalInteg(FluidNF,vfes,pfes);
   res->SetVDotGF(vdot_gf);
   if(windkessel_bc || resistance_bc) 
   {
      res->AddCVBC_Op(cvbc_op,&A3d_diag,&Abc_diag,&Atotal_diag);
   }

   // immersed fsi to residual operator
   if (IFEM)
   {
      res->SetIFEM(*im_solid,f_fsi);
   }

   if (direct_solver) res->UsingDirectSolver();
   if (monolithic_prec) res->UsingMonolithicPrec();

   res->SetTimers(&sw_fluidtan,&sw_fluidres);

   //////////////////////////////////////////////////////////////////////////////////
   // setup newton solver
   #ifdef MFEM_USE_MPI
   newton_solver = new FluidNewtonSolver(res,block_trueOffsets,vfes->GetComm());
   #else
   newton_solver = new FluidNewtonSolver(res,block_trueOffsets);
   #endif

   // weak bcs
   newton_solver->WeakBCs(&vel_nbcs,L_L_form,ft,&accel_terms,f_form,fb);

   // solution vectors
   newton_solver->SetSolutionVectors(yn,ydotn,ynp1,ydotnp1,yg,ydotg);
   if (total_integ) newton_solver->SetTotalInteg();

   // tell newton solver about the time stepping method
   newton_solver->UseGenAlpha(alpha_f,alpha_m,gamma);

   if (direct_solver) newton_solver->UsingDirectSolver();

   // cv bcs
   if (windkessel_bc) newton_solver->AddWindkesselBC(&windkessel_bcs);

   // immersed fsi
   if (IFEM)
   {
      newton_solver->SetIFEM(*im_solid);
   }

   // flag to not abort when newton solver does not converge
   // this is useful for immersed fsi with thin structures and
   // nonzero scaling of tau_m to reduce inaccuracies from poor
   // approximation of pressure across the interface
   if (ns_no_error_nonconvergence)
   {
      newton_solver->SetNoErrorNonConvergence();
   }

   // timing
   sw_tan.Clear();
   sw_res.Clear();
   sw_linsolve.Clear();
   newton_solver->SetTimers(&sw_tan,&sw_res,&sw_linsolve, &sw_initiators, &sw_correctors, &sw_weakbcs);

   // options for newton solver
   newton_solver->SetMinIter(minIt_n);
   newton_solver->iterative_mode = true;

   if (no_prec)
   {
      if (zero_mean)
      {
         #ifndef MFEM_USE_MPI
         solver_wrap = new BlockOrthoSolver(block_trueOffsets);
         #else
         solver_wrap = new BlockOrthoSolver(block_trueOffsets,vfes->GetComm());
         #endif
         solver_wrap->SetSolver(*solver); 

         newton_solver->SetSolver(*solver_wrap);
      }
      else
      {
         newton_solver->SetSolver(*solver);
      }
   }
   else
   {
      if (zero_mean)
      {
         // zero pressure mean after each newton iteration if necessary
         #ifndef MFEM_USE_MPI
         solver_wrap = new BlockOrthoSolver(block_trueOffsets);
         #else
         solver_wrap = new BlockOrthoSolver(block_trueOffsets,vfes->GetComm());
         #endif
         solver_wrap->SetSolver(*jacobian_prec);

         static_cast<ExtendedFGMRESSolver *>(solver)->SetPreconditioner(*solver_wrap);
      }
      else
      {
         static_cast<ExtendedFGMRESSolver *>(solver)->SetPreconditioner(*jacobian_prec);
      }

      newton_solver->SetSolver(*solver);
   }  

   newton_solver->SetOperator(*res);
   newton_solver->SetPrintLevel(pl_n);
   newton_solver->SetRelTol(rtol_n);
   newton_solver->SetAbsTol(atol_n);
   newton_solver->SetMaxIter(maxIt_n);
   newton_solver->PrintNewton(print_newton);
   
   /// INITIAL CONDITIONS ///
   // If the initial condition was set, it has to be aligned with dependent
   // Vectors and GridFunctions
   Vector &vn = yn.GetBlock(0);
   Vector &pn = yn.GetBlock(1);
   #ifdef MFEM_USE_MPI
   vn_gf.ParallelProject(vn);
   pn_gf.ParallelProject(pn);
   #else
   vn_gf.GetTrueDofs(vn);
   pn_gf.GetTrueDofs(pn);
   #endif

   sw_setup.Stop();
}

void IncompressibleFluidSolver::Step(double &time, double dt, int current_step,
                        bool provisional)
{
   if (verbose && my_rank == 0)
   {
      if (current_step == 1)
      {
         mfem::out << "---------------------------------------" << std::endl;
         mfem::out << "             TIME STEPPING             " << std::endl;
         mfem::out << "---------------------------------------" << std::endl;
      }

      if (restarting && current_step == restart_step + 1)
      {
         mfem::out << "--------------------------------------------------------------" << std::endl;
         mfem::out << "RESTARTING SIMULATION FROM:" << std::endl;
         mfem::out << "    TIME = " << time - dt << std::endl;
         mfem::out << "    STEP = " << restart_step << std::endl;
      }

      if (current_step != 1)
      {
         mfem::out << "--------------------------------------------------------------" << std::endl;
      }
      mfem::out << std::setw(6) << "Step " << std::setw(6) << current_step << std::setw(8) << " Time " << std::setw(8) << " dt " << "\n";
      
      mfem::out << std::setw(15) << "           " << std::setw(3) << std::setprecision(2)
                << std::scientific << time << "  " << std::setw(4) << dt << "\n";
   }

   sw_step.Start();

   // Update parameters //
   UpdateParameters(time,dt,current_step); 

   // Predictor Step //
   Predictors(dt);

   // Set Essential BCs //
   SetEssentialBCs(time); // time = t_n+1

   // Newton Iterations //
   sw_solve.Start();
   newton_solver->NewtonSolve(dt,time);
   sw_solve.Stop();
   iter_solve = newton_solver->GetNumIterations();
   res_solve = newton_solver->GetFinalNorm();
   rel_res_solve = newton_solver->GetFinalRelNorm();

   // Get residuals
   static_cast<FluidNewtonSolver *>(newton_solver)->GetResiduals(norm,rel_norm,norm_mom,rel_norm_mom,norm_mass,rel_norm_mass);

   // Update variables //
   UpdateVariables(time,dt);

   // Output to screen and save data //

   // Write restarts
   if (write_restarts && (current_step % restart_save_step == 0) && !provisional)
   {
      writeRestarts(current_step,dt,time);
   }

   sw_step.Stop();
}

bool IncompressibleFluidSolver::Step(double &time, double dt, int current_step)
{
   if (verbose && my_rank == 0)
   {
      if (current_step == 1)
      {
         mfem::out << "---------------------------------------" << std::endl;
         mfem::out << "             TIME STEPPING             " << std::endl;
         mfem::out << "---------------------------------------" << std::endl;
      }

      if (restarting && current_step == restart_step + 1)
      {
         mfem::out << "--------------------------------------------------------------" << std::endl;
         mfem::out << "RESTARTING SIMULATION FROM:" << std::endl;
         mfem::out << "    TIME = " << time - dt << std::endl;
         mfem::out << "    STEP = " << restart_step << std::endl;
      }

      if (current_step != 1)
      {
         mfem::out << "--------------------------------------------------------------" << std::endl;
      }
      mfem::out << std::setw(6) << "Step " << std::setw(6) << current_step << std::setw(8) << " Time " << std::setw(8) << " dt " << "\n";
      
      mfem::out << std::setw(15) << "           " << std::setw(3) << std::setprecision(2)
                << std::scientific << time << "  " << std::setw(4) << dt << "\n";
   }

   sw_step.Start();

   // Update parameters //
   UpdateParameters(time,dt,current_step); 

   // Predictor Step //
   Predictors(dt);

   // Set Essential BCs //
   SetEssentialBCs(time); // time = t_n+1

   // Newton Iterations //
   sw_solve.Start();
   newton_solver->NewtonSolve(dt,time);
   sw_solve.Stop();
   if (newton_solver->GetConverged() == false)
   {
      return false;
   }
   iter_solve = newton_solver->GetNumIterations();
   res_solve = newton_solver->GetFinalNorm();
   rel_res_solve = newton_solver->GetFinalRelNorm();

   // Get residuals
   static_cast<FluidNewtonSolver *>(newton_solver)->GetResiduals(norm,rel_norm,norm_mom,rel_norm_mom,norm_mass,rel_norm_mass);

   sw_step.Stop();

   return true;
}

void IncompressibleFluidSolver::UpdateParameters(double time, double dt, int current_step_)
{
   cur_step = current_step_;
   curr_dt = dt;

   if (IFEM)
   {
      im_solid->SetTimeAndDt(time,dt);
   }
}

void IncompressibleFluidSolver::Predictors(double dt)
{
   sw_predictors.Start();
   ynp1 = yn;
   ydotnp1.Set((gamma - 1.0)/gamma,ydotn);

   if (zero_mean)
   {
      // remove mean from pressure
      Vector &pnp1 = ynp1.GetBlock(1);
      Orthogonalize(pnp1);
   }

   // immersed solid
   if (IFEM)
   {
      Vector &vnp1 = ynp1.GetBlock(0);
      Vector &vdotnp1 = ydotnp1.GetBlock(0);

      im_solid->Predictors(vnp1,vdotnp1);
   }
   sw_predictors.Stop();
}

void IncompressibleFluidSolver::SetEssentialBCs(const double &time)
{
   sw_essbcs.Start();
   Vector &vnp1 = ynp1.GetBlock(0);
   Vector &pnp1 = ynp1.GetBlock(1);

   // velocity
   vn_next_gf.SetFromTrueDofs(vnp1);
   for (auto &vel_dbc : vel_dbcs)
   {
      if(vel_dbc.all_comp == true)
      {
         vel_dbc.coeff->SetTime(time);
         vn_next_gf.ProjectBdrCoefficient(*vel_dbc.coeff,vel_dbc.attr);
      }
      else
      {
         for (int d = 0; d < pmesh->Dimension(); d++)
         {
            if (vel_dbc.coeff_array[d] == nullptr)
            {
               continue;
            }
            else
            {
               vel_dbc.coeff_array[d]->SetTime(time);
            }
         }
         vn_next_gf.ProjectBdrCoefficient(vel_dbc.coeff_array, vel_dbc.attr);
      }
   }
   #ifdef MFEM_USE_MPI
   vn_next_gf.ParallelProject(vnp1);
   #else
   vn_next_gf.GetTrueDofs(vnp1);
   #endif
   // zero boundary nodes of a prescribed velocity bc
   for (auto &vel_dbc : vel_dbcs)
   {
      if(vel_dbc.zero_bdr)
      {
         const int N = vel_dbc.zero_bdr_tdofs.Size();
         const auto idx = vel_dbc.zero_bdr_tdofs.Read();
         auto Y_RW = vnp1.ReadWrite();
         mfem::forall(N, [=] MFEM_HOST_DEVICE (int i) { Y_RW[idx[i]] = 0.0; });
      }
      else
      {
         continue;
      }
   }

   // pressure
   pn_next_gf.SetFromTrueDofs(pnp1);
   for (auto &pres_dbc : pres_dbcs)
   {
      pres_dbc.coeff->SetTime(time);
      pn_next_gf.ProjectBdrCoefficient(*pres_dbc.coeff,pres_dbc.attr);
   }
   // individual pressure node
   if (pressure_node_condition)
   {
      if (p_local_node_num > -1)
      {
         if (ref_p_coeff)
         {
            Vector coord(pmesh->Dimension());
            coord = 0.0; // not used
            double ref_p_val = ref_p_coeff(coord,time);
            pn_next_gf(p_local_node_num) = ref_p_val;
         }
         else
         {
            pn_next_gf(p_local_node_num) = ref_p;
         }
      }
   }
   #ifdef MFEM_USE_MPI
   pn_next_gf.ParallelProject(pnp1);
   #else
   pn_next_gf.GetTrueDofs(pnp1);   
   #endif

   // Zero out essential dofs in vdot and pdot for Gen-alpha scheme
   Vector &vdot_np1 = ydotnp1.GetBlock(0);
   Vector &pdot_np1 = ydotnp1.GetBlock(1);

   // velocity
   const int N = vel_ess_tdofs.Size();
   const auto idx = vel_ess_tdofs.Read();
   auto Y_RW = vdot_np1.ReadWrite();
   mfem::forall(N, [=] MFEM_HOST_DEVICE (int i) { Y_RW[idx[i]] = 0.0; });

   // pressure
   const int Np = pres_ess_tdofs.Size();
   const auto idxp = pres_ess_tdofs.Read();
   auto Y_RWp = pdot_np1.ReadWrite();
   mfem::forall(Np, [=] MFEM_HOST_DEVICE (int i) { Y_RWp[idxp[i]] = 0.0; });

   if (IFEM)
   {      
      // NOTE: fluid should also have no slip condition satisfied on the
      // fixed immersed boundaries

      im_solid->SetEssentialBCs(time);
   }
   sw_essbcs.Stop();
}

void IncompressibleFluidSolver::UpdateVariables(double &time, double &dt)
{
   sw_updates.Start();
   if (windkessel_bc) UpdateWindkesselParameters();
   
   UpdateTimestepHistory(time-dt,dt);

   if (IFEM) im_solid->Update();
   sw_updates.Stop();
}

void IncompressibleFluidSolver::UpdateTimestepHistory(const double &time, double &dt)
{
   Vector &vn = yn.GetBlock(0);
   Vector &vnp1 = ynp1.GetBlock(0);

   Vector &pn = yn.GetBlock(1);
   Vector &pnp1 = ynp1.GetBlock(1);

   if (zero_mean)
   {
      pn_next_gf.SetFromTrueDofs(pnp1);
      MeanZero(pn_next_gf);
      #ifdef MFEM_USE_MPI
      pn_next_gf.ParallelProject(pnp1);
      #else
      pn_next_gf.GetTrueDofs(pnp1);
      #endif
   }

   Vector &vdotn = ydotn.GetBlock(0);
   Vector &vdotnp1 = ydotnp1.GetBlock(0);

   Vector &pdotn = ydotn.GetBlock(1);
   Vector &pdotnp1 = ydotnp1.GetBlock(1);

   // velocity
   vn = vnp1;
   vn_gf.SetFromTrueDofs(vn);

   vdotn = vdotnp1;
   vdotn_gf.SetFromTrueDofs(vdotn);

   // pressure
   pn = pnp1;
   pn_gf.SetFromTrueDofs(pn);

   pdotn = pdotnp1;
   pdotn_gf.SetFromTrueDofs(pdotn);
}

void IncompressibleFluidSolver::UpdateWindkesselParameters()
{
   // compute current integ term with new converged velocity, then update 
   // prev time, prev flow rate, and integ_hist_terms and pass integ_hist_terms to each 
   // windkessel boundary

   for (auto &wbc : windkessel_bcs)
   {
      // compute hist terms
      double temp = 0.0;

      // compute wbc.q_np1 //
      // q_np1 = (q_npaf - qn)/af + qn;
      wbc.q_np1 = (wbc.q_npalphaf - wbc.q_n) / alpha_f + wbc.q_n;

      // midpoint rule approx
      double t_mid = 0.5*(wbc.t_np1 + wbc.t_n);
      temp = wbc.dt * (exp(t_mid/(wbc.Rd*wbc.C)) / wbc.C ) * (((t_mid - wbc.t_n)/(wbc.t_np1 - wbc.t_n)) * wbc.q_np1 
                           + ((t_mid - wbc.t_np1)/(wbc.t_n - wbc.t_np1)) * wbc.q_n);
      wbc.hist_terms += temp;

      // set previous variables accordingly
      // time
      wbc.t_nm1 = wbc.t_n;
      wbc.t_n = wbc.t_np1;

      // flow rates
      wbc.q_nm1 = wbc.q_n;
      wbc.q_n = wbc.q_np1;
   }
}

void IncompressibleFluidSolver::PrintInfo()
{
   #ifdef MFEM_USE_MPI
   int fes_size0 = vfes->GlobalTrueVSize();
   int fes_size1 = pfes->GlobalTrueVSize();
   int numelements = pmesh->GetNE();
   MPI_Allreduce(MPI_IN_PLACE, &numelements, 1, MPI_INT, MPI_SUM, vfes->GetComm());
   #else
   int fes_size0 = vfes->GetTrueVSize();
   int fes_size1 = pfes->GetTrueVSize();
   int numelements = pmesh->GetNE();
   #endif

   if (my_rank == 0)
   {
      mfem::out <<   "IFS version: " << IFS_VERSION << std::endl
                  << "MFEM version: " << MFEM_VERSION << std::endl
                  << "Number of fluid elements: " << numelements << std::endl
                  << "Velocity #DOFs: " << fes_size0 << std::endl
                  << "Pressure #DOFs: " << fes_size1 << std::endl
                  << "Total #DOFs: " << fes_size1 + fes_size0 << std::endl
                  << "Velocity order = " << orderv << std::endl
                  << "Pressure order = " << orderp << std::endl;
   }
}

void IncompressibleFluidSolver::RegisterTimers()
{
   if (timers_registered) return;
   timer_reg.Clear();

   // Base rows
   timer_reg.Add("setup",     "SETUP",     &sw_setup,     "timestep", "", 0, false);
   timer_reg.Add("timestep",  "TIMESTEP",  &sw_step,      "timestep", "", 0, false);

   if (extended_timing)
   {
      timer_reg.Add("predictors","PREDICTORS",&sw_predictors,"timestep", "", 1);
      timer_reg.Add("essbcs",    "ESSBCS",    &sw_essbcs,    "timestep", "", 1);
   }

   // Level-1 row
   timer_reg.Add("newton",    "NEWTONSOLVE", &sw_solve,   "timestep", "", 1);

   if (extended_timing)
   {
      timer_reg.Add("initiators","INITIATORS",&sw_initiators,"timestep", "newton", 2);
      timer_reg.Add("weakbcs",   "WEAKBCS",   &sw_weakbcs,   "timestep", "newton", 2);
   }

   // RESIDUAL and children
   timer_reg.Add("residual",  "RESIDUAL",    &sw_res,     "timestep", "newton", 2);
   if (IFEM)
   {
      timer_reg.Add("fluid_res", "FLUIDRES", &sw_fluidres, "timestep", "residual", 3);
      timer_reg.Add("solid_res", "SOLIDRES", &sw_solidres, "timestep", "residual", 3);
      timer_reg.Add("fsi_res",   "FSIRES",   &sw_fsires,   "timestep", "residual", 3);
   }

   // TANGENT and children
   timer_reg.Add("tangent",   "TANGENT",     &sw_tan,     "timestep", "newton", 2);
   if (IFEM)
   {
      timer_reg.Add("fluid_tan", "FLUIDTAN", &sw_fluidtan, "timestep", "tangent", 3);
      timer_reg.Add("solid_tan", "SOLIDTAN", &sw_solidtan, "timestep", "tangent", 3);
      timer_reg.Add("fsi_tan",   "FSITAN",   &sw_fsitan,   "timestep", "tangent", 3);
   }

   // FSICMAT under NEWTONSOLVE (same level as TANGENT/RESIDUAL/LINSOLVE)
   if (IFEM)
   {
      timer_reg.Add("fsi_cmat", "FSICMAT", &sw_fsicmat, "timestep", "newton", 2);
   }

   // LINSOLVE (kept at same level under NEWTONSOLVE)
   timer_reg.Add("linsolve",  "LINSOLVE",    &sw_linsolve,"timestep", "newton", 2);

   if (extended_timing)
   {
      timer_reg.Add("correctors","CORRECTORS",&sw_correctors,"timestep", "newton", 2);
   }

   if (extended_timing) timer_reg.Add("updates","UPDATES", &sw_updates,"timestep","",1);
   if (write_restarts) timer_reg.Add("restarts","RESTARTS",&sw_restarts,"timestep","",1);
   timer_reg.Add("io","I/O",&sw_io,"timestep","",0,false);

   timers_registered = true;
}

void IncompressibleFluidSolver::PrintTimingData()
{
   RegisterTimers();

   // Collect raw times in seconds
   std::vector<double> times = timer_reg.CollectSeconds();

#ifdef MFEM_USE_MPI
   // Average across ranks
   std::vector<double> sum(times.size());
   MPI_Allreduce(times.data(), sum.data(),
                 static_cast<int>(times.size()),
                 MPI_DOUBLE, MPI_SUM, vfes->GetComm());
   int num_procs = Mpi::WorldSize();
   for (size_t i = 0; i < times.size(); ++i) times[i] = sum[i] / static_cast<double>(num_procs);
#endif

   const auto& entries = timer_reg.Entries();
   auto time_of = [&](const std::string& key) {
      for (size_t i=0; i<entries.size(); ++i)
         if (entries[i].key == key) return times[i];
      return 0.0;
   };

   if (my_rank == 0)
   {
      mfem::out << "\n";
      mfem::out << "--------------------------TIMER INFO--------------------------\n";
      mfem::out << std::setprecision(2);
      mfem::out << std::setw(5) << "FORMAT = " << std::setw(26) << "DAY-HR:MIN:SEC.MILISECONDs\n\n";

      std::ostringstream oss;
      oss << std::fixed << std::setprecision(2);

      // Total label field width (indent + label + padding) before the time column.
      // Keep this constant so the time column is vertically aligned.
      const int label_field_total = 20;  // tune as needed
      const double base_t = time_of("timestep");

      auto print_row = [&](const std::string& indent, const std::string& label,
                           const char* time_str, const std::string& pct_str,
                           bool show_pct)
      {
         // pad between label end and time column so that time_str starts at fixed column
         int used = static_cast<int>(indent.size()) + static_cast<int>(label.size());
         int pad = std::max(1, label_field_total - used);
         mfem::out << indent
                   << std::left << label
                   << std::string(pad, ' ')
                   << std::right << std::setw(13) << time_str;
         if (show_pct)
         {
            mfem::out << " (" << pct_str << "% of TIMESTEP)";
         }
         mfem::out << "\n";
      };

      for (size_t i = 0; i < entries.size(); ++i)
      {
         const auto& e = entries[i];
         if (!e.enabled) continue;

         // Adjust the reported time for certain timers to avoid double
         // counting any time accumulated prior to a restart. For example,
         // for the SETUP timer we subtract e.restart_time here while the
         // raw times[] vector still contains the full cumulative time.
         double this_time = times[i];
         if (e.key == "setup")
         {
            this_time -= e.restart_time;
            if (this_time < 0.0) { this_time = 0.0; }
         }

         char buf[32];
         ConvertTimeToFormattedTime(this_time, buf, sizeof(buf));

         std::string indent;
         if (e.indent > 0) indent.assign(3 * e.indent, ' ');

         double denom = e.percent_of.empty() ? base_t : time_of(e.percent_of);
         double pct = (denom > 0.0) ? (this_time / denom * 100.0) : 0.0;
         oss.str("");
         oss.clear();
         oss << std::setw(5) << pct;

         print_row(indent, e.label, buf, oss.str(), e.show_percent);
      }

      mfem::out << std::setprecision(8);
   }
}

void IncompressibleFluidSolver::ConvertTimeToFormattedTime(double time_in_seconds, char* time_str, size_t size)
{
   // Decompose using milliseconds to avoid floating rounding issues.
   long long total_ms = static_cast<long long>(std::llround(time_in_seconds * 1000.0));
   long long ms = total_ms % 1000LL;
   long long total_seconds = total_ms / 1000LL;

   long long days = total_seconds / 86400LL;
   int hours = static_cast<int>((total_seconds % 86400LL) / 3600LL);
   int minutes = static_cast<int>((total_seconds % 3600LL) / 60LL);
   int seconds = static_cast<int>(total_seconds % 60LL);

   // Format as: day-HR:MIN:SEC.MILISECONDs
   // Example: 0-00:00:00.000
   std::snprintf(time_str, size, "%lld-%02d:%02d:%02d.%03lld",
                 days, hours, minutes, seconds, ms);
}

///////////////////////////
/// BOUNDARY CONDITIONS ///
///////////////////////////

void IncompressibleFluidSolver::AddVelDirichletBC(VectorCoefficient *coeff, Array<int> &attr, bool zero_bdr)
{
   vel_dbcs.emplace_back(attr, coeff);

   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Velocity Dirichlet BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }

   int dim = pmesh->Dimension();

   for (int i = 0; i < attr.Size(); ++i)
   {
      for (int j = 0; j < pmesh->Dimension(); j++)
      {
         // check each component on this bdr attribute
         MFEM_ASSERT((vel_ess_attr[i*dim + j] && attr[i]) == 0,
                  "Duplicate boundary definition detected for component " + std::to_string(j) + ".");
      }

      // only if no duplicate for all components then add into vel_ess_attr
      for (int j = 0; j < pmesh->Dimension(); j++)
      {
         if (attr[i] == 1)
         {
            vel_ess_attr[i*dim + j] = 1;
         }
      }
   }

   // add essential dofs to list
   Array<int> temp_dof;
   vfes->GetEssentialTrueDofs(attr,temp_dof);
   vel_ess_tdofs.Append(temp_dof,temp_dof.Size());

   if (zero_bdr)
   {
      vel_dbcs.back().zero_bdr = true;

      // extract all the true dofs besides this bc
      Array<int> bdr_tdofs;
      Array<int> not_attr;
      not_attr.SetSize(attr.Size());
      not_attr = 1;
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            not_attr[i] = 0;
         }
      }
      vfes->GetEssentialTrueDofs(not_attr,bdr_tdofs);
      bdr_tdofs.Sort();

      // Find intersection of bdr_tdofs and temp_dof (essential dofs for this bc)
      Array<int> intersect_tdofs;
      for (int i=0; i < temp_dof.Size(); i++)
      {
         int idx = bdr_tdofs.FindSorted(temp_dof[i]);
         if (idx >= 0)
         {
            intersect_tdofs.Append(temp_dof[i]);
         }
      }
      vel_dbcs.back().zero_bdr_tdofs = intersect_tdofs;
   }
}

void IncompressibleFluidSolver::AddVelDirichletBC(VecFuncT *f, Array<int> &attr, bool zero_bdr)
{
   AddVelDirichletBC(new VectorFunctionCoefficient(pmesh->Dimension(), f), attr, zero_bdr);
}

void IncompressibleFluidSolver::AddVelDirichletBC(Coefficient *coeff, Array<int> &attr, int component_, bool zero_bdr)
{
   if (verbose && my_rank == 0)
   {
      std::string comp_str;
      if (component_ == 1) comp_str = "X";
      else if (component_ == 2) comp_str = "Y";
      else if (component_ == 3) comp_str = "Z";
      else MFEM_ABORT("Invalid velocity component (Choose 1 = X, 2 = Y, or 3 = Z)!");
      mfem::out << "Adding Velocity Dirichlet BC to component " << comp_str << " and attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }

   int component = component_ - 1;

   // form coefficient array, and put null if unset and coefficient if set inside temp
   std::unique_ptr<Coefficient*[]> temp(new Coefficient*[pmesh->Dimension()]); // pointer to an array of pointers to Coefficient objects of size dimension
   for(int j= 0; j < pmesh->Dimension(); j++)
   {
      if (component == j)
      {
         temp[j] = coeff;
      }
      else
      {
         temp[j] = nullptr;
      }
   }

   int dim = pmesh->Dimension();
   vel_dbcs.emplace_back(attr, temp.release(), component);
   for (int i = 0; i < attr.Size(); ++i)
   {
      // check each component on this bdr attribute
      MFEM_ASSERT((vel_ess_attr[i*dim + component] && attr[i]) == 0,
               "Duplicate boundary definition detected for component " + std::to_string(component_) + ".");

      if (attr[i] == 1)
      {
         vel_ess_attr[i*dim+component] = 1;
      }
   }

   // add essential dofs to list
   Array<int> temp_dof;
   vfes->GetEssentialTrueDofs(attr,temp_dof,component);
   vel_ess_tdofs.Append(temp_dof,temp_dof.Size());

   if (zero_bdr)
   {
      vel_dbcs.back().zero_bdr = true;

      // extract all the true dofs besides this bc
      Array<int> bdr_tdofs;
      Array<int> not_attr;
      not_attr.SetSize(attr.Size());
      not_attr = 1;
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            not_attr[i] = 0;
         }
      }
      vfes->GetEssentialTrueDofs(not_attr,bdr_tdofs,component);
      bdr_tdofs.Sort();

      // Find intersection of bdr_tdofs and temp_dof (essential dofs for this bc)
      Array<int> intersect_tdofs;
      for (int i=0; i < temp_dof.Size(); i++)
      {
         int idx = bdr_tdofs.FindSorted(temp_dof[i]);
         if (idx >= 0)
         {
            intersect_tdofs.Append(temp_dof[i]);
         }
      }
      vel_dbcs.back().zero_bdr_tdofs = intersect_tdofs;
   }
}

void IncompressibleFluidSolver::AddVelDirichletBC(ScalarFuncT *f, Array<int> &attr, int component_, bool zero_bdr)
{
   AddVelDirichletBC(new FunctionCoefficient(f), attr, component_, zero_bdr);
}

void IncompressibleFluidSolver::AddPresDirichletBC(Coefficient *coeff, Array<int> &attr)
{
   MFEM_ABORT("Error: Dont use these! Use AddPressureBC instead!");
   pres_dbcs.emplace_back(attr, coeff);

   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Pressure Dirichlet BC to attributes ";
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
      MFEM_ASSERT((pres_ess_attr[i] && attr[i]) == 0,
                  "Duplicate boundary definition deteceted.");
      if (attr[i] == 1)
      {
         pres_ess_attr[i] = 1;
      }
   }

   zero_mean = false;

   // add essential dofs to list
   Array<int> temp_dof;
   pfes->GetEssentialTrueDofs(attr,temp_dof);
   pres_ess_tdofs.Append(temp_dof,temp_dof.Size());
}

void IncompressibleFluidSolver::AddPresDirichletBC(ScalarFuncT *f, Array<int> &attr)
{
   AddPresDirichletBC(new FunctionCoefficient(f), attr);
}

void IncompressibleFluidSolver::AddPressureNodeCondition(const double &p, Vector &coords, double tol)
{
   if (verbose && my_rank == 0)
   {
      mfem::out << "Setting Nodal Pressure Value to ";
      mfem::out << p << " at Node Position " << "(" << coords(0) << "," << coords(1);
      if (pmesh->Dimension() == 3) mfem::out << "," << coords(2);
      mfem::out << ")";
      mfem::out << std::endl;
   }

   zero_mean = false;
   pressure_node_condition = true;
   ref_p = p;
   ref_p_coeff = nullptr;

   // find the node number with position provided by coords
   DenseMatrix point_mat(pmesh->Dimension());
   point_mat.SetCol(0,coords);
   Array<int> elem_ids;
   Array<IntegrationPoint> ips;
   int npts_found = pmesh->FindPoints(point_mat,elem_ids,ips);
   MFEM_ASSERT(npts_found == 1, "Error: node to set pressure value was not found!");
   int elem_num = elem_ids[0];
   p_local_node_num = -1;
   int dim = pmesh->Dimension();
   Vector curr_ip(dim);
   if (elem_num > -1)
   {
      const FiniteElement *pfe = pfes->GetFE(elem_num);
      ElementTransformation *eltrans = pfes->GetElementTransformation(elem_num);
      Array<int> vdofs;
      auto doftrans = pfes->GetElementVDofs(elem_num,vdofs);
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

         // come nodes with coords
         diff_x = abs(curr_ip(0) - coords(0));
         diff_y = abs(curr_ip(1) - coords(1));
         if (dim == 3) diff_z = abs(curr_ip(2) - coords(2));

         if (diff_x <= tol && diff_y <= tol)
         {
            if (dim == 3)
            {
               if (diff_z <= tol)
               {
                  p_local_node_num = vdofs[j];
                  break;
               }
            }
            else
            {
               p_local_node_num = vdofs[j];
               break;
            }
         }
      }

      MFEM_VERIFY(p_local_node_num > -1, "Error: node not found!");

      // add found node to the list
      pres_ess_tdofs.Append(p_local_node_num);
   }
}

void IncompressibleFluidSolver::AddPressureNodeCondition(ScalarFuncT &p_coeff, Vector &coords, double tol)
{
   if (verbose && my_rank == 0)
   {
      mfem::out << "Setting Nodal Pressure Value";
      mfem::out << " at Node Position " << "(" << coords(0) << "," << coords(1);
      if (pmesh->Dimension() == 3) mfem::out << "," << coords(2);
      else mfem::out << ")";
      mfem::out << std::endl;
   }

   zero_mean = false;
   pressure_node_condition = true;
   ref_p_coeff = &p_coeff;

   // find the node number with position provided by coords
   DenseMatrix point_mat(pmesh->Dimension());
   point_mat.SetCol(0,coords);
   Array<int> elem_ids;
   Array<IntegrationPoint> ips;
   pmesh->FindPoints(point_mat,elem_ids,ips);
   int elem_num = elem_ids[0];
   p_local_node_num = -1;
   int dim = pmesh->Dimension();
   Vector curr_ip(dim);
   if (elem_num > -1)
   {
      const FiniteElement *pfe = pfes->GetFE(elem_num);
      ElementTransformation *eltrans = pfes->GetElementTransformation(elem_num);
      Array<int> vdofs;
      auto doftrans = pfes->GetElementVDofs(elem_num,vdofs);
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

         // come nodes with coords
         double diff_x = abs(curr_ip(0) - coords(0));
         double diff_y = abs(curr_ip(1) - coords(1));
         if (dim == 3) double diff_z = abs(curr_ip(2) - coords(2));

         if (diff_x <= tol && diff_y <= tol)
         {
            if (dim == 3)
            {
               if (diff_z <= tol)
               {
                  p_local_node_num = vdofs[j];
                  break;
               }
            }
            else
            {
               p_local_node_num = vdofs[j];
               break;
            }
         }
      }

      MFEM_VERIFY(p_local_node_num > -1, "Error: node not found!");

      // add found node to the list
      pres_ess_tdofs.Append(p_local_node_num);
   }
}

void IncompressibleFluidSolver::AddPeriodicBC(Array<int> &attr)
{
   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Periodic BC to attributes ";
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
      MFEM_ASSERT((periodic_attr[i] && attr[i]) == 0,
                  "Duplicate boundary definition deteceted.");
      if (attr[i] == 1)
      {
         periodic_attr[i] = 1;
      }
   }
}

void IncompressibleFluidSolver::AddVelNeumannBC(VectorCoefficient *coeff, Array<int> &attr)
{
   vel_nbcs.emplace_back(attr, coeff);

   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Traction BC to attributes ";
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
}

void IncompressibleFluidSolver::AddVelNeumannBC(VecFuncT *f, Array<int> &attr)
{
   AddVelNeumannBC(new VectorFunctionCoefficient(pmesh->Dimension(), f), attr);
}

void IncompressibleFluidSolver::AddVelNeumannBC(MatrixCoefficient *coeff, Array<int> &attr)
{
   vel_nbcs.emplace_back(attr, coeff);

   if (verbose)
   {
      mfem::out << "Adding Traction BC to attributes ";
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
}

void IncompressibleFluidSolver::AddVelNeumannBC(MatFuncT *f, Array<int> &attr)
{
   AddVelNeumannBC(new MatrixFunctionCoefficient(pmesh->Dimension(), f), attr);
}

void IncompressibleFluidSolver::AddVelNeumannBC(Coefficient *coeff, Array<int> &attr, bool pressure)
{
   vel_nbcs.emplace_back(attr, coeff,  pressure);

   if (verbose && my_rank == 0)
   {
      if (pressure)
      {
         mfem::out << "Adding Pressure BC (Normal Traction) to attributes ";
      }
      else
      {
         mfem::out << "Adding Traction BC to attributes ";
      }
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
}

void IncompressibleFluidSolver::AddVelNeumannBC(ScalarFuncT *f, Array<int> &attr, bool pressure)
{
   AddVelNeumannBC(new FunctionCoefficient(f), attr, pressure);
}

void IncompressibleFluidSolver::AddOutflowBC(Array<int> &attr)
{  
   AddTractionFreeBC(attr);
}

void IncompressibleFluidSolver::AddNoSlipBC(Array<int> &attr)
{
   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding No Slip BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }

   // AddVelDirichletBC(neumann_vel,attr);
   vel_dbcs.emplace_back(attr, new VectorFunctionCoefficient(pmesh->Dimension(),&no_slip));
   int dim = pmesh->Dimension();

   for (int i = 0; i < attr.Size(); ++i)
   {
      for (int j = 0; j < pmesh->Dimension(); j++)
      {
         // check each component on this bdr attribute
         MFEM_ASSERT((vel_ess_attr[i*dim + j] && attr[i]) == 0,
                  "Duplicate boundary definition detected for component " + std::to_string(j) + ".");
      }

      // only if no duplicate for all components then add into vel_ess_attr
      for (int j = 0; j < pmesh->Dimension(); j++)
      {
         if (attr[i] == 1)
         {
            vel_ess_attr[i*dim + j] = 1;
         }
      }
   }

   // add essential dofs to list
   Array<int> temp_dof;
   vfes->GetEssentialTrueDofs(attr,temp_dof);
   vel_ess_tdofs.Append(temp_dof,temp_dof.Size());
}

void IncompressibleFluidSolver::AddBodyForce(VectorCoefficient *coeff, Array<int> &attr)
{
   accel_terms.emplace_back(attr, coeff);

   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Body Force to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }
}

void IncompressibleFluidSolver::AddBodyForce(VecFuncT *f, Array<int> &attr)
{
   AddBodyForce(new VectorFunctionCoefficient(pmesh->Dimension(), f), attr);
}

void IncompressibleFluidSolver::UseBackflowStabilization(Array<int> &attr, double beta_)
{
   backflow_stab = true;

   // duplicate bc check
   for (auto &bf_bc : backflow_stabilization_bcs)
   {
      for (int i = 0; i < attr.Size(); ++i)
      {
         MFEM_ASSERT((bf_bc.attr[i] && attr[i]) == 0,
                     "Duplicate boundary definition detected.");
      }
   }

   // if no duplicate, continue to add bc  
   backflow_stabilization_bcs.emplace_back(attr,beta_);


   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Backflow Stabilization to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }
}

// resistance bc
#ifdef MFEM_USE_MPI
void IncompressibleFluidSolver::AddResistanceBC(Array<int> &attr, double R_, double p0_, ParFiniteElementSpace *gfes)
#else
void IncompressibleFluidSolver::AddResistanceBC(Array<int> &attr, double R_, double p0_, FiniteElementSpace *gfes)
#endif
{
   // duplicate bc check
   for (auto &rbc : resistance_bcs)
   {
      for (int i = 0; i < attr.Size(); ++i)
      {
         MFEM_ASSERT((rbc.attr[i] && attr[i]) == 0,
                     "Duplicate boundary definition detected.");
      }
   }

   // if no duplicate, continue to add bc
   resistance_bcs.emplace_back(attr, R_, p0_, gfes ? gfes : vfes);

   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Resistance BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }   

   // set resistance bc flag to true
   resistance_bc = true;
   zero_mean = false;
}

// windkessel bc
#ifdef MFEM_USE_MPI
void IncompressibleFluidSolver::AddWindkesselBC(Array<int> &attr, double R_, double Rd_, double C_, Coefficient *p_0_, std::function<double(double)> p_d_, double init_flow_rate_, ParFiniteElementSpace *gfes)
#else
void IncompressibleFluidSolver::AddWindkesselBC(Array<int> &attr, double R_, double Rd_, double C_, Coefficient *p_0_, std::function<double(double)> p_d_, double init_flow_rate_, FiniteElementSpace *gfes)
#endif
{
   // duplicate bc check
   for (auto &wbc : windkessel_bcs)
   {
      for (int i = 0; i < attr.Size(); ++i)
      {
         MFEM_ASSERT((wbc.attr[i] && attr[i]) == 0,
                     "Duplicate boundary definition detected.");
      }
   }

   // if no duplicate, continue to add bc
   windkessel_bcs.emplace_back(attr, R_, Rd_, C_, p_0_, p_d_, init_flow_rate_, gfes ? gfes : vfes);

   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Windkessel BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }   

   // set windkessel bc flag to true
   windkessel_bc = true;
   zero_mean = false;
}

#ifdef MFEM_USE_MPI
void IncompressibleFluidSolver::AddWindkesselBC(Array<int> &attr, double R_, double Rd_, double C_, double p_0_, std::function<double(double)> p_d_, double init_flow_rate_, ParFiniteElementSpace *gfes)
#else
void IncompressibleFluidSolver::AddWindkesselBC(Array<int> &attr, double R_, double Rd_, double C_, double p_0_, std::function<double(double)> p_d_, double init_flow_rate_, FiniteElementSpace *gfes)
#endif
{
   // duplicate bc check
   for (auto &wbc : windkessel_bcs)
   {
      for (int i = 0; i < attr.Size(); ++i)
      {
         MFEM_ASSERT((wbc.attr[i] && attr[i]) == 0,
                     "Duplicate boundary definition detected.");
      }
   }

   // if no duplicate, continue to add bc  
   windkessel_bcs.emplace_back(attr, R_, Rd_, C_, p_0_, p_d_, init_flow_rate_, gfes ? gfes : vfes);  

   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Windkessel BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }   

   // set windkessel bc flag to true
   windkessel_bc = true;
   zero_mean = false;
}

#ifdef MFEM_USE_MPI
void IncompressibleFluidSolver::AddWindkesselBC(Array<int> &attr, double R_, double Rd_, double C_, double p_0_, double p_d_, double init_flow_rate_, ParFiniteElementSpace *gfes)
#else
void IncompressibleFluidSolver::AddWindkesselBC(Array<int> &attr, double R_, double Rd_, double C_, double p_0_, double p_d_, double init_flow_rate_, FiniteElementSpace *gfes)
#endif
{
   // duplicate bc check
   for (auto &wbc : windkessel_bcs)
   {
      for (int i = 0; i < attr.Size(); ++i)
      {
         MFEM_ASSERT((wbc.attr[i] && attr[i]) == 0,
                     "Duplicate boundary definition detected.");
      }
   }

   // if no duplicate, continue to add bc  
   windkessel_bcs.emplace_back(attr, R_, Rd_, C_, p_0_, p_d_, init_flow_rate_, gfes ? gfes : vfes);   

   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Windkessel BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }   

   // set windkessel bc flag to true
   windkessel_bc = true;
   zero_mean = false;
}

#ifdef MFEM_USE_MPI
void IncompressibleFluidSolver::AddMappedInflowBC_FromFile(int inflow_type_, Mesh *mesh, std::string filename_, bool time_depend, Array<int> &attr, int approx_type, bool periodic_data, int num_fourier_modes_)
{
   mapped_inflow = new MappedInflowCoefficient(pmesh->Dimension(),true,inflow_type_,orderv,pmesh,mesh,attr,time_depend);
   mapped_inflow->FlowRateFileOptions(filename_, approx_type, periodic_data, num_fourier_modes_);
   
   vel_dbcs.emplace_back(attr, mapped_inflow);

   if (verbose && pmesh->GetMyRank() == 0)
   {
      mfem::out << "Adding Mapped Inflow Velocity Dirichlet BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }

   // for (int i = 0; i < attr.Size(); ++i)
   // {
   //    MFEM_ASSERT((vel_ess_attr[i] && attr[i]) == 0,
   //                "Duplicate boundary definition detected.");
   //    if (attr[i] == 1)
   //    {
   //       vel_ess_attr[i] = 1;
   //    }
   // }

   int dim = pmesh->Dimension();

   for (int i = 0; i < attr.Size(); ++i)
   {
      for (int j = 0; j < pmesh->Dimension(); j++)
      {
         // check each component on this bdr attribute
         MFEM_ASSERT((vel_ess_attr[i*dim + j] && attr[i]) == 0,
                  "Duplicate boundary definition detected for component " + std::to_string(j) + ".");
      }

      // only if no duplicate for all components then add into vel_ess_attr
      for (int j = 0; j < pmesh->Dimension(); j++)
      {
         if (attr[i] == 1)
         {
            vel_ess_attr[i*dim + j] = 1;
         }
      }
   }

   // add essential dofs to list
   Array<int> temp_dof;
   vfes->GetEssentialTrueDofs(attr,temp_dof);
   vel_ess_tdofs.Append(temp_dof,temp_dof.Size());
}

void IncompressibleFluidSolver::AddMappedInflowBC(int inflow_type_, Mesh *mesh, std::function<double (double)> flow_rate, bool time_depend, Array<int> &attr)
{
   mapped_inflow = new MappedInflowCoefficient(pmesh->Dimension(),flow_rate,inflow_type_,orderv,pmesh,mesh,attr,time_depend);
   
   vel_dbcs.emplace_back(attr, mapped_inflow);

   if (verbose && pmesh->GetMyRank() == 0)
   {
      mfem::out << "Adding Mapped Inflow Velocity Dirichlet BC to attributes ";
      for (int i = 0; i < attr.Size(); ++i)
      {
         if (attr[i] == 1)
         {
            mfem::out << i << " ";
         }
      }
      mfem::out << std::endl;
   }

   // for (int i = 0; i < attr.Size(); ++i)
   // {
   //    MFEM_ASSERT((vel_ess_attr[i] && attr[i]) == 0,
   //                "Duplicate boundary definition detected.");
   //    if (attr[i] == 1)
   //    {
   //       vel_ess_attr[i] = 1;
   //    }
   // }

   int dim = pmesh->Dimension();

   for (int i = 0; i < attr.Size(); ++i)
   {
      for (int j = 0; j < pmesh->Dimension(); j++)
      {
         // check each component on this bdr attribute
         MFEM_ASSERT((vel_ess_attr[i*dim + j] && attr[i]) == 0,
                  "Duplicate boundary definition detected for component " + std::to_string(j) + ".");
      }

      // only if no duplicate for all components then add into vel_ess_attr
      for (int j = 0; j < pmesh->Dimension(); j++)
      {
         if (attr[i] == 1)
         {
            vel_ess_attr[i*dim + j] = 1;
         }
      }
   }

   // add essential dofs to list
   Array<int> temp_dof;
   vfes->GetEssentialTrueDofs(attr,temp_dof);
   vel_ess_tdofs.Append(temp_dof,temp_dof.Size());
}
#endif

#ifdef MFEM_USE_MPI
void IncompressibleFluidSolver::AddImmersedSolid(ParMesh *solid_pmesh, int orderd, Coefficient *init_rhos, MaterialType *mat_, int form_type_, bool visc_, double scale_fac_)
#else
void IncompressibleFluidSolver::AddImmersedSolid(Mesh *solid_pmesh, int orderd, Coefficient *init_rhos, MaterialType *mat_, int form_type_, bool visc_, double scale_fac_)
#endif
{
   if (verbose && my_rank == 0)
   {
      mfem::out << "Adding Immersed Solid";
      // TO DO: name here
      mfem::out << std::endl;
   }

   IFEM = true;
   im_solid = new ImmersedSolid(solid_pmesh,orderd,init_rhos,fluid_density,mu_f.constant,&accel_terms,vfes,pmesh,mat_,orderingtype,form_type_,visc_);

   if (f_fsi.Size() == 0)
   {
      // initialize
      f_fsi.SetSize(vfes->GetTrueVSize());
      f_fsi = 0.0;

      scale_fac = scale_fac_;
      if(scale_fac != 1.0)
      {
         s_gf.SetSpace(pfes);
         s_gf = 0.0;

         im_solid->StabsScaling(&s_gf);
      }      
   }

   MFEM_ASSERT(dimensional_flag, "Error: fsi only supported with dimensional form of equations!");
}

void IncompressibleFluidSolver::AddImmersedSolid(ImmersedSolid *im_solid_, double scale_fac_)
{
   if (verbose)
   {
      mfem::out << "Adding Immersed Solid";
      // TO DO: name here
      mfem::out << std::endl;
   }

   IFEM = true;
   
   im_solid = im_solid_;

   if (f_fsi.Size() == 0)
   {
      // initialize
      f_fsi.SetSize(vfes->GetTrueVSize());
      f_fsi = 0.0;

      scale_fac = scale_fac_;
      if(scale_fac != 1.0)
      {
         s_gf.SetSpace(pfes);
         s_gf = 0.0;

         im_solid->StabsScaling(&s_gf);
      }      
   }

   MFEM_ASSERT(dimensional_flag, "Error: fsi only supported with dimensional form of equations!");
}

///////////////////////
/// POST-PROCESSING ///
///////////////////////

void IncompressibleFluidSolver::ComputeCurl3D(GridFunction &u, GridFunction &cu)
{
   FiniteElementSpace *fes = u.FESpace();

   // AccumulateAndCountZones.
   Array<int> zones_per_vdof;
   zones_per_vdof.SetSize(fes->GetVSize());
   zones_per_vdof = 0;

   cu = 0.0;

   // Local interpolation.
   int elndofs;
   Array<int> vdofs;
   Vector vals;
   Vector loc_data;
   int vdim = fes->GetVDim();
   DenseMatrix grad_hat;
   DenseMatrix dshape;
   DenseMatrix grad;
   Vector curl;

   for (int e = 0; e < fes->GetNE(); ++e)
   {
      fes->GetElementVDofs(e, vdofs);
      u.GetSubVector(vdofs, loc_data);
      vals.SetSize(vdofs.Size());
      ElementTransformation *tr = fes->GetElementTransformation(e);
      const FiniteElement *el = fes->GetFE(e);
      elndofs = el->GetDof();
      int dim = el->GetDim();
      dshape.SetSize(elndofs, dim);

      for (int dof = 0; dof < elndofs; ++dof)
      {
         // Project.
         const IntegrationPoint &ip = el->GetNodes().IntPoint(dof);
         tr->SetIntPoint(&ip);

         // Eval and GetVectorGradientHat.
         el->CalcDShape(tr->GetIntPoint(), dshape);
         grad_hat.SetSize(vdim, dim);
         DenseMatrix loc_data_mat(loc_data.GetData(), elndofs, vdim);
         MultAtB(loc_data_mat, dshape, grad_hat);

         const DenseMatrix &Jinv = tr->InverseJacobian();
         grad.SetSize(grad_hat.Height(), Jinv.Width());
         Mult(grad_hat, Jinv, grad);

         curl.SetSize(3);
         curl(0) = grad(2, 1) - grad(1, 2);
         curl(1) = grad(0, 2) - grad(2, 0);
         curl(2) = grad(1, 0) - grad(0, 1);

         for (int j = 0; j < curl.Size(); ++j)
         {
            vals(elndofs * j + dof) = curl(j);
         }
      }

      // Accumulate values in all dofs, count the zones.
      for (int j = 0; j < vdofs.Size(); j++)
      {
         int ldof = vdofs[j];
         cu(ldof) += vals[j];
         zones_per_vdof[ldof]++;
      }
   }

   #ifdef MFEM_USE_MPI
   // Communication
   // Count the zones globally.
   GroupCommunicator &gcomm = dynamic_cast<ParFiniteElementSpace *>(u.FESpace())->GroupComm();
   gcomm.Reduce<int>(zones_per_vdof, GroupCommunicator::Sum);
   gcomm.Bcast(zones_per_vdof);

   // Accumulate for all vdofs.
   gcomm.Reduce<double>(cu.GetData(), GroupCommunicator::Sum);
   gcomm.Bcast<double>(cu.GetData());
   #endif

   // Compute means.
   for (int i = 0; i < cu.Size(); i++)
   {
      const int nz = zones_per_vdof[i];
      if (nz)
      {
         cu(i) /= nz;
      }
   }
}

void IncompressibleFluidSolver::ComputeCurl2D(GridFunction &u,
                                 GridFunction &cu,
                                 bool assume_scalar)
{
   FiniteElementSpace *fes = u.FESpace();

   // AccumulateAndCountZones.
   Array<int> zones_per_vdof;
   zones_per_vdof.SetSize(fes->GetVSize());
   zones_per_vdof = 0;

   cu = 0.0;

   // Local interpolation.
   int elndofs;
   Array<int> vdofs;
   Vector vals;
   Vector loc_data;
   int vdim = fes->GetVDim();
   DenseMatrix grad_hat;
   DenseMatrix dshape;
   DenseMatrix grad;
   Vector curl;

   for (int e = 0; e < fes->GetNE(); ++e)
   {
      fes->GetElementVDofs(e, vdofs);
      u.GetSubVector(vdofs, loc_data);
      vals.SetSize(vdofs.Size());
      ElementTransformation *tr = fes->GetElementTransformation(e);
      const FiniteElement *el = fes->GetFE(e);
      elndofs = el->GetDof();
      int dim = el->GetDim();
      dshape.SetSize(elndofs, dim);

      for (int dof = 0; dof < elndofs; ++dof)
      {
         // Project.
         const IntegrationPoint &ip = el->GetNodes().IntPoint(dof);
         tr->SetIntPoint(&ip);

         // Eval and GetVectorGradientHat.
         el->CalcDShape(tr->GetIntPoint(), dshape);
         grad_hat.SetSize(vdim, dim);
         DenseMatrix loc_data_mat(loc_data.GetData(), elndofs, vdim);
         MultAtB(loc_data_mat, dshape, grad_hat);

         const DenseMatrix &Jinv = tr->InverseJacobian();
         grad.SetSize(grad_hat.Height(), Jinv.Width());
         Mult(grad_hat, Jinv, grad);

         if (assume_scalar)
         {
            curl.SetSize(2);
            curl(0) = grad(0, 1);
            curl(1) = -grad(0, 0);
         }
         else
         {
            curl.SetSize(2);
            curl(0) = grad(1, 0) - grad(0, 1);
            curl(1) = 0.0;
         }

         for (int j = 0; j < curl.Size(); ++j)
         {
            vals(elndofs * j + dof) = curl(j);
         }
      }

      // Accumulate values in all dofs, count the zones.
      for (int j = 0; j < vdofs.Size(); j++)
      {
         int ldof = vdofs[j];
         cu(ldof) += vals[j];
         zones_per_vdof[ldof]++;
      }
   }

   #ifdef MFEM_USE_MPI
   // Communication.
   // Count the zones globally.
   GroupCommunicator &gcomm = dynamic_cast<ParFiniteElementSpace *>(u.FESpace())->GroupComm();
   gcomm.Reduce<int>(zones_per_vdof, GroupCommunicator::Sum);
   gcomm.Bcast(zones_per_vdof);

   // Accumulate for all vdofs.
   gcomm.Reduce<double>(cu.GetData(), GroupCommunicator::Sum);
   gcomm.Bcast<double>(cu.GetData());
   #endif

   // Compute means.
   for (int i = 0; i < cu.Size(); i++)
   {
      const int nz = zones_per_vdof[i];
      if (nz)
      {
         cu(i) /= nz;
      }
   }
}

double IncompressibleFluidSolver::ComputeCFL(const GridFunction &u, double dt)
{
   Mesh *pmesh_u = u.FESpace()->GetMesh();
   const FiniteElementSpace *fes = u.FESpace();
   int vdim = fes->GetVDim();

   Vector ux, uy, uz;
   Vector ur, us, ut;
   double cflx = 0.0;
   double cfly = 0.0;
   double cflz = 0.0;
   double cflm = 0.0;
   double cflmax = 0.0;

   for (int e = 0; e < fes->GetNE(); ++e)
   {
      const FiniteElement *fe = fes->GetFE(e);
      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                               fe->GetOrder());
      ElementTransformation *tr = fes->GetElementTransformation(e);

      u.GetValues(e, ir, ux, 1);
      ur.SetSize(ux.Size());
      u.GetValues(e, ir, uy, 2);
      us.SetSize(uy.Size());
      if (vdim == 3)
      {
         u.GetValues(e, ir, uz, 3);
         ut.SetSize(uz.Size());
      }

      double hmin = pmesh_u->GetElementSize(e, 1) /
                    (double) fes->GetElementOrder(0);

      for (int i = 0; i < ir.GetNPoints(); ++i)
      {
         // const IntegrationPoint &ip = ir.IntPoint(i);
         // tr->SetIntPoint(&ip);
         // const DenseMatrix &invJ = tr->InverseJacobian();
         // const double detJinv = 1.0 / tr->Jacobian().Det();

         // if (vdim == 2)
         // {
         //    ur(i) = (ux(i) * invJ(0, 0) + uy(i) * invJ(1, 0)) * detJinv;
         //    us(i) = (ux(i) * invJ(0, 1) + uy(i) * invJ(1, 1)) * detJinv;
         // }
         // else if (vdim == 3)
         // {
         //    ur(i) = (ux(i) * invJ(0, 0) + uy(i) * invJ(1, 0)
         //             + uz(i) * invJ(2, 0))
         //            * detJinv;
         //    us(i) = (ux(i) * invJ(0, 1) + uy(i) * invJ(1, 1)
         //             + uz(i) * invJ(2, 1))
         //            * detJinv;
         //    ut(i) = (ux(i) * invJ(0, 2) + uy(i) * invJ(1, 2)
         //             + uz(i) * invJ(2, 2))
         //            * detJinv;
         // }

         cflx = fabs(dt * ux(i) / hmin);
         cfly = fabs(dt * uy(i) / hmin);
         if (vdim == 3)
         {
            cflz = fabs(dt * uz(i) / hmin);
         }
         cflm = cflx + cfly + cflz;
         cflmax = fmax(cflmax, cflm);
      }
   }
   
   #ifdef MFEM_USE_MPI
   MPI_Allreduce(MPI_IN_PLACE, &cflmax, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, vfes->GetComm());
   #endif

   return cflmax;
}

double IncompressibleFluidSolver::ComputeIntegralDivV(const GridFunction &v_gf)
{
   const FiniteElement *fe;
   ElementTransformation *Tr;
   DenseMatrix gradv;
   int intorder;
   int dim = v_gf.FESpace()->GetMesh()->SpaceDimension();
   double div_v;
   double div_int = 0.0;
   gradv.SetSize(dim);
   gradv = 0.0;

   double p = 1.0;

   for (int i = 0; i < v_gf.FESpace()->GetNE(); i++)
   {
      fe = v_gf.FESpace()->GetFE(i);
      Tr = v_gf.FESpace()->GetElementTransformation(i);
      intorder = 2*fe->GetOrder() + 3; // <--------
      const IntegrationRule *ir;
      ir = &(IntRules.Get(fe->GetGeomType(), intorder));
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir->IntPoint(j);
         Tr->SetIntPoint(&ip);
         v_gf.GetVectorGradient(*Tr,gradv);

         // compute div v
         div_v = gradv.Trace();

         // integral
         div_int += ip.weight * Tr->Weight() * div_v;
      }
   }

   #ifdef MFEM_USE_MPI
   MPI_Allreduce(MPI_IN_PLACE, &div_int, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   #endif

   return div_int;
}

// method for adaptive time stepping based on CFL_Crit
double IncompressibleFluidSolver::estimate_convective_cfl_dt()
{
   int vdim = vfes->GetVDim();

   Vector vx, vy, vz;
   real_t dt_x = mfem::infinity();
   real_t dt_y = mfem::infinity();
   real_t dt_z = mfem::infinity();
   real_t dt_min_local = mfem::infinity();
   real_t dt_min = mfem::infinity();
   real_t u_eps = 1e-10;

   for (int e = 0; e < vfes->GetNE(); ++e)
   {
      const FiniteElement *fe = vfes->GetFE(e);
      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                               fe->GetOrder());
      ElementTransformation *tr = vfes->GetElementTransformation(e);

      vn_gf.GetValues(e, ir, vx, 1);
      vn_gf.GetValues(e, ir, vy, 2);
      if (vdim == 3)
      {
         vn_gf.GetValues(e, ir, vz, 3);
      }

      real_t hmin = pmesh->GetElementSize(e, 1) /
                    (real_t) orderv;

      for (int i = 0; i < ir.GetNPoints(); ++i)
      {
         dt_x = CFL_crit * hmin / (fabs(vx(i)) + u_eps);
         dt_y = CFL_crit * hmin / (fabs(vy(i)) + u_eps);
         if (vdim == 3)
         {
            dt_z = CFL_crit * hmin / (fabs(vz(i)) + u_eps);
         }

         dt_min_local = std::min({dt_x,dt_y,dt_z});
         if (dt_min_local < dt_min) 
         {
            dt_min = dt_min_local;
         }      
      }
   }

   #ifdef MFEM_USE_MPI
   real_t dtmin_global = 0.0;
   MPI_Allreduce(&dt_min,
                 &dtmin_global,
                 1,
                 MPITypeMap<real_t>::mpi_type,
                 MPI_MIN,
                 vfes->GetComm());
   #elif
   real_t dtmin_global = dt_min;
   #endif

   return dtmin_global;
}

void IncompressibleFluidSolver::ComputeWSS(GridFunction &u, GridFunction &p, int avg_type)
{
   /*
   // Wall shear stress = force that the fluid exerts on the wall

   t = sigma * n
   WSS_fluid = t - (t*n)n --> shear force that fluid experiences (force that wall exerts on fluid)

   WSS = -WSS_fluid
   n_wall = -n; (n is outward normal)
   t_wall = sigma * n_wall
   WSS = t_wall - (t_wall * n_wall) n_wall


   // sigma*n computed using the fluid velocity/pressure = force that wall exerts on fluid (or the force that the fluid experiences)
   // to get the force that fluid exerts on the wall (or force the wall experiences) multiply by minus sign
   // Think of pressure driven channel, flow left to right. Fluid flows and experiences an opposing friction (force wall exerts on fluid) and
   // computing wall shear stress in this case will be negative!
   // Note: here we just average values across element boundaries (cheap and decent results, SPR procedure would be better)
   */

   if (wss_gf.FESpace() == nullptr) wss_gf.SetSpace(vfes);

   if (avg_type == ifs_solver::WSSAvgType::Standard)
   {
      bool dimensional_eq = false;
      if (dimensional_flag == ifs_solver::dim_flag::Dimensional) dimensional_eq = true;

      // zero gf
      wss_gf = 0.0;

      // AccumulateAndCountZones
      Array<int> zones_per_vdof;
      zones_per_vdof.SetSize(vfes->GetVSize());
      zones_per_vdof = 0;

      // loop over wall elements and compute
      const FiniteElement *ufe;
      FaceElementTransformations *ftr;
      ElementTransformation *eltrans;
      FiniteElementSpace *ufes = u.FESpace();

      // Local interpolation
      int dim = pmesh->Dimension();
      Vector nor;
      bool on_surf = false;
      DenseMatrix eye(dim), gradu(dim), sigma(dim);
      Vector wss_elvec, tract(dim);
      Array<int> vdofs;

      // loop over boundary elements
      for (int e = 0; e < pmesh -> GetNBE(); e++)
      {         
         // only compute if on surface of interest (TO DO)
         on_surf = true;
         if(on_surf) 
         {         
            eltrans = pmesh->GetBdrElementTransformation(e);
            ufe = ufes->GetBE(e);
            int dim = ufe->GetDim() + 1; // bdr elment for cg are 1 dim less than actual
            auto doftrans = ufes -> GetBdrElementVDofs (e, vdofs);

            nor.SetSize(vdim);
            nor = 0.0;

            int nd = ufe->GetDof();
            wss_elvec.SetSize(nd*dim);

            // integration rule = nodes of element (velocity)
            // i.e. int point 0 = node 0
            const IntegrationRule *ir = &ufe->GetNodes();

            for (int j = 0; j < ir->GetNPoints(); j++)
            {
               // setup current integration point
               const IntegrationPoint &ip = ir->IntPoint(j);
               eltrans->SetIntPoint(&ip);

               // compute n // compute outward normal
               if (dim > 1)
               {
                  CalcOrtho(eltrans->Jacobian(), nor);
               }
               else
               {
                  nor[0] = 1.0;
               }
               // normalize nor, see example 18 in general it is not a unit normal
               double normag = 0;
               for (int ii = 0; ii < nor.Size(); ii++)
               {
                  normag += nor(ii)*nor(ii);
               }
               normag = sqrt(normag);
               nor *= 1/normag;

               // to get WSS the wall experiences
               nor *= -1.0;

               // compute stress
               // Compute -p*I
               eye.Diag(1.0,vdim);
               eye *= -1.0 * p.GetValue(*eltrans,ip);

               // -p*I + mu*(gradu + gradu^T)
               u.GetVectorGradient(*eltrans,gradu);
               sigma = gradu;
               gradu.Transpose();
               sigma += gradu; // gradu^T
               sigma *= dimensional_eq ? kin_vis * fluid_density : kin_vis; // mu = rho * nu
               sigma += eye;

               // trac = sigma * n
               sigma.Mult(nor,tract);

               // shear trac = trac - (trac * n) n
               double normal_tract = tract*nor;
               // store at the node
               wss_elvec(j) = tract(0) - normal_tract * nor(0);
               wss_elvec(j + nd) = tract(1) - normal_tract * nor(1);
               if (dim == 3) wss_elvec(j + 2*nd) = tract(2) - normal_tract * nor(2);
            }

            // Accumulate values in all dofs, count the zones.
            for (int j = 0; j < vdofs.Size(); j++)
            {
               int ldof = vdofs[j];
               wss_gf(ldof) += wss_elvec[j];
               zones_per_vdof[ldof]++;
            }
         }
      }

      #ifdef MFEM_USE_MPI
      // Communication
      // Count the zones globally.
      GroupCommunicator &gcomm = wss_gf.ParFESpace()->GroupComm();
      gcomm.Reduce<int>(zones_per_vdof, GroupCommunicator::Sum);
      gcomm.Bcast(zones_per_vdof);

      // Accumulate for all vdofs.
      gcomm.Reduce<double>(wss_gf.GetData(), GroupCommunicator::Sum);
      gcomm.Bcast<double>(wss_gf.GetData());
      #endif

      // Compute means
      for (int i = 0; i < wss_gf.Size(); i++)
      {
         const int nz = zones_per_vdof[i];
         if (nz)
         {
            wss_gf(i) /= nz;
         }
      }
   }
   else if (avg_type == ifs_solver::WSSAvgType::LumpedL2)
   {
      MFEM_ASSERT(equal_order,"Known issue with mixed order spaces");
      bool dimensional_eq = false;
      if (dimensional_flag == ifs_solver::dim_flag::Dimensional) dimensional_eq = true;

      // zero gf
      wss_gf = 0.0;

      // AccumulateAndCountZones
      Array<int> zones_per_vdof;
      zones_per_vdof.SetSize(vfes->GetVSize());
      zones_per_vdof = 0;

      Vector shape_gf(vfes->GetVSize());
      shape_gf = 0.0;

      // loop over wall elements and compute
      const FiniteElement *ufe;
      FaceElementTransformations *ftr;
      ElementTransformation *eltrans;
      FiniteElementSpace *ufes = u.FESpace();

      // Local interpolation
      int dim = pmesh->Dimension();
      Vector nor;
      bool on_surf = false;
      DenseMatrix eye(dim), gradu(dim), sigma(dim);
      Vector wss_elvec, shape_elvec, tract(dim);
      Array<int> vdofs;

      // loop over boundary elements
      for (int e = 0; e < pmesh -> GetNBE(); e++)
      {         
         // only compute if on surface of interest (TO DO)
         on_surf = true;
         if(on_surf) 
         {         
            eltrans = pmesh->GetBdrElementTransformation(e);
            ufe = ufes->GetBE(e);
            int dim = ufe->GetDim() + 1; // bdr elment for cg are 1 dim less than actual
            auto doftrans = ufes -> GetBdrElementVDofs (e, vdofs);

            nor.SetSize(vdim);
            nor = 0.0;

            int nd = ufe->GetDof();
            Vector shape(nd);
            wss_elvec.SetSize(nd*dim);
            shape_elvec.SetSize(nd*dim);

            wss_elvec = 0.0;
            shape_elvec = 0.0;

            // integration rule
            int integ_order = 2*ufe->GetOrder() + eltrans->OrderW();
            const IntegrationRule *ir = &IntRules.Get(eltrans->GetGeometryType(), integ_order);

            for (int j = 0; j < ir->GetNPoints(); j++)
            {
               // setup current integration point
               const IntegrationPoint &ip = ir->IntPoint(j);
               eltrans->SetIntPoint(&ip);

               double w = ip.weight * eltrans->Weight();

               // compute shape functions
               ufe->CalcShape(ip,shape);

               // compute n // compute outward normal
               if (dim > 1)
               {
                  CalcOrtho(eltrans->Jacobian(), nor);
               }
               else
               {
                  nor[0] = 1.0;
               }
               // normalize nor, see example 18 in general it is not a unit normal
               double normag = 0;
               for (int ii = 0; ii < nor.Size(); ii++)
               {
                  normag += nor(ii)*nor(ii);
               }
               normag = sqrt(normag);
               nor *= 1/normag;

               // to get WSS the wall experiences
               nor *= -1.0;

               // compute stress
               // Compute -p*I
               eye.Diag(1.0,vdim);
               eye *= -1.0 * p.GetValue(*eltrans,ip);

               // -p*I + mu*(gradu + gradu^T)
               u.GetVectorGradient(*eltrans,gradu);
               sigma = gradu;
               gradu.Transpose();
               sigma += gradu; // gradu^T
               sigma *= dimensional_eq ? kin_vis * fluid_density : kin_vis; // mu = rho * nu
               sigma += eye;

               // trac = sigma * n
               sigma.Mult(nor,tract);

               // shear trac = trac - (trac * n) n
               double normal_tract = tract*nor;
               // store at the node
               for (int a = 0; a < nd; a++)
               {
                  wss_elvec(a) += w * shape(a) * (tract(0) - normal_tract * nor(0));
                  wss_elvec(a + nd) += w * shape(a) * (tract(1) - normal_tract * nor(1));
                  if (dim == 3) wss_elvec(a + 2*nd) += w * shape(a) * (tract(2) - normal_tract * nor(2));

                  shape_elvec(a) += w * shape(a);
                  shape_elvec(a + nd) += w * shape(a);
                  if (dim == 3) shape_elvec(a + 2*nd) += w * shape(a);
               }               
            }
            
            // Accumulate values in all dofs, count the zones.
            for (int j = 0; j < vdofs.Size(); j++)
            {
               int ldof = vdofs[j];
               wss_gf(ldof) += wss_elvec[j];
               shape_gf(ldof) += shape_elvec[j];
            }
         }
      }

      #ifdef MFEM_USE_MPI
      // Communication
      // Count the zones globally.
      GroupCommunicator &gcomm = wss_gf.ParFESpace()->GroupComm();
      gcomm.Reduce<int>(zones_per_vdof, GroupCommunicator::Sum);
      gcomm.Bcast(zones_per_vdof);

      // Accumulate for all vdofs.
      gcomm.Reduce<double>(wss_gf.GetData(), GroupCommunicator::Sum);
      gcomm.Bcast<double>(wss_gf.GetData());
      #endif

      // Compute dof values
      for (int i = 0; i < wss_gf.Size(); i++)
      {
         const double shape_val = shape_gf[i];
         if (abs(shape_val) > 0.0)
         {
            wss_gf(i) /= shape_val;
         }
      }
   }
   else
   {
      MFEM_ABORT("Error: undefined averaging type for WSS");
   }   
}

// computes mass (or volume flow rate) thru a surface
double IncompressibleFluidSolver::ComputeMassFlowThruSurf(GridFunction &u, Array<int> &attr_surf, bool inflow_flag, double fluid_density, bool volume_flow_rate)
{
   // \int_{Gamma} \rho u dot n dGamma (by default)
   // or 
   // \int_{Gamma} u dot n dGamma

   // or in ale case
   
   // \int_{Gamma} \rho (u-uhat) dot n dGamma (by default)
   // or 
   // \int_{Gamma} (u-uhat) dot n dGamma


   const FiniteElement *ufe;
   FaceElementTransformations *ftr;
   ElementTransformation *eltrans;
   FiniteElementSpace *ufes = u.FESpace();

   // Local interpolation
   int dim = pmesh->Dimension();
   Vector vel_val, nor, uhat_val;
   int vdim = ufes->GetVDim();
   vel_val.SetSize(vdim);
   vel_val = 0.0;
   uhat_val.SetSize(vdim);
   uhat_val = 0.0;
   double integ = 0.0;
   bool on_surf = false;

   // loop over boundary elements
   for (int e = 0; e < pmesh -> GetNBE(); e++)
   {
      on_surf = false;

      // check boundary attribute is on surface of interest
      int bdr_attr = pmesh->GetBdrAttribute(e);

      // Note: bdr attributes have same number as specified in GMSH, i.e bdr 1 = bdr attribute 1 in MFEM
      // just when you set the attribute arrays is the bdr attribute need -1;
      for (int ii = 0; ii < attr_surf.Size(); ii++)
      {
         if(attr_surf[ii] == 1 && bdr_attr == ii+1) // +1 bc of above comment
         {
            on_surf = true;
         }
      }
      
      // only compute if on surface of interest
      if(on_surf) 
      {
         eltrans = pmesh->GetBdrElementTransformation(e);
         ufe = ufes->GetBE(e);
         int dim = ufe->GetDim() + 1; // bdr elment for cg are 1 dim less than actual

         nor.SetSize(vdim);
         nor = 0.0;

         // integration rule
         // exact GLL integration for 2*ctx.order integrand (computing on finer mesh)
         int order;
         order = ufe->GetOrder() + 1;
         const IntegrationRule *ir = &IntRules.Get(eltrans->GetGeometryType(), order);

         for (int j = 0; j < ir->GetNPoints(); j++)
         {
            // setup current integration point
            const IntegrationPoint &ip = ir->IntPoint(j);
            eltrans->SetIntPoint(&ip);

            // compute u
            u.GetVectorValue(*eltrans,ip,vel_val);
            if (ale)
            {
               MFEM_ABORT("TO DO");
               // // compute uhat
               // mesh_velocity.GetVectorValue(*eltrans,ip,uhat_val);
               // vel_val -= uhat_val;
            }

            // compute n // compute outward normal
            if (dim > 1)
            {
               CalcOrtho(eltrans->Jacobian(), nor);
            }
            else
            {
               nor[0] = 1.0;
            }
            // normalize nor, see example 18 in general it is not a unit normal
            double normag = 0;
            for (int ii = 0; ii < nor.Size(); ii++)
            {
               normag += nor(ii)*nor(ii);
            }
            normag = sqrt(normag);
            nor *= 1/normag;
            if (inflow_flag) { nor *= -1.0; } // to get mass flow in
            
            // rho * u dot n // compute integral
            double val = (vel_val * nor) * ip.weight * eltrans->Weight();
            if (!volume_flow_rate) { val *= fluid_density;}

            integ += val;
         }
      }
   }

   #ifdef MFEM_USE_MPI
   MPI_Allreduce(MPI_IN_PLACE, &integ, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   #endif
   
   return integ;
}

// computes mass (or volume flow rate) thru a surface
double IncompressibleFluidSolver::ComputeAvgVeloOnSurf(GridFunction &u, Array<int> &attr_surf, Vector &avg_velo)
{
   const FiniteElement *ufe;
   FaceElementTransformations *ftr;
   ElementTransformation *eltrans;
   FiniteElementSpace *ufes = u.FESpace();

   // Local interpolation
   Vector vel_val;
   int vdim = ufes->GetVDim();
   vel_val.SetSize(vdim);
   vel_val = 0.0;

   double integ_ux = 0.0;
   double integ_uy = 0.0;
   double integ_uz = 0.0;
   double integ_area = 0.0;
   bool on_surf = false;

   // loop over boundary elements
   for (int e = 0; e < pmesh -> GetNBE(); e++)
   {
      on_surf = false;

      // check boundary attribute is on surface of interest
      int bdr_attr = pmesh->GetBdrAttribute(e);

      // Note: bdr attributes have same number as specified in GMSH, i.e bdr 1 = bdr attribute 1 in MFEM
      // just when you set the attribute arrays is the bdr attribute need -1;
      for (int ii = 0; ii < attr_surf.Size(); ii++)
      {
         if(attr_surf[ii] == 1 && bdr_attr == ii+1) // +1 bc of above comment
         {
            on_surf = true;
         }
      }
      
      // only compute if on surface of interest
      if(on_surf) 
      {
         eltrans = pmesh->GetBdrElementTransformation(e);
         ufe = ufes->GetBE(e);
         int dim = ufe->GetDim() + 1; // bdr elment for cg are 1 dim less than actual

         // integration rule
         // exact GLL integration for 2*ctx.order integrand (computing on finer mesh)
         int order;
         order = ufe->GetOrder() + 1;
         const IntegrationRule *ir = &IntRules.Get(eltrans->GetGeometryType(), order);

         for (int j = 0; j < ir->GetNPoints(); j++)
         {
            // setup current integration point
            const IntegrationPoint &ip = ir->IntPoint(j);
            eltrans->SetIntPoint(&ip);

            // compute u
            u.GetVectorValue(*eltrans,ip,vel_val);

            integ_ux += ip.weight * eltrans->Weight() * vel_val(0);
            integ_uy += ip.weight * eltrans->Weight() * vel_val(1);
            if (dim == 3) integ_uz += ip.weight * eltrans->Weight() * vel_val(2);

            integ_area += ip.weight * eltrans->Weight();
         }
      }
   }

   #ifdef MFEM_USE_MPI
   MPI_Allreduce(MPI_IN_PLACE, &integ_ux, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &integ_uy, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   if(pmesh->Dimension() == 3) MPI_Allreduce(MPI_IN_PLACE, &integ_uz, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &integ_area, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   #endif

   MFEM_ASSERT(integ_area > 0.0, "Error: zero area for surface");
   avg_velo(0) = integ_ux / integ_area;
   avg_velo(1) = integ_uy / integ_area;
   if (pmesh->Dimension() == 3) avg_velo(2) = integ_uz / integ_area;
   
   return avg_velo.Norml2();
}

// computes average pressure on surface
double IncompressibleFluidSolver::ComputeAvgPressureOnSurf(GridFunction &p_gf, Array<int> &attr_surf, double &area)
{
   // \int_{Gamma} p dGamma / \int{Gamma} dGamma

   const FiniteElement *pfe;
   ElementTransformation *eltrans;
   FiniteElementSpace *pfes = p_gf.FESpace();

   // Local interpolation
   double integ_p = 0.0;
   double integ_area = 0.0;
   bool on_surf = false;

   // loop over boundary elements
   for (int e = 0; e < pmesh -> GetNBE(); e++)
   {
      on_surf = false;

      // check boundary attribute is on surface of interest
      int bdr_attr = pmesh->GetBdrAttribute(e);

      // Note: bdr attributes have same number as specified in GMSH, i.e bdr 1 = bdr attribute 1 in MFEM
      // just when you set the attribute arrays is the bdr attribute need -1;
      for (int ii = 0; ii < attr_surf.Size(); ii++)
      {
         if(attr_surf[ii] == 1 && bdr_attr == ii+1) // +1 bc of above comment
         {
            on_surf = true;
         }
      }
      
      // only compute if on surface of interest
      if(on_surf) 
      {
         eltrans = pmesh->GetBdrElementTransformation(e);
         pfe = pfes->GetBE(e);
         int dim = pfe->GetDim() + 1; // bdr elment for cg are 1 dim less than actual

         // integration rule
         // exact GLL integration for ctx.order integrand (computing on finer mesh)
         int order;
         order = pfe->GetOrder() + 1;
         const IntegrationRule *ir = &IntRules.Get(eltrans->GetGeometryType(), order);

         for (int j = 0; j < ir->GetNPoints(); j++)
         {
            // setup current integration point
            const IntegrationPoint &ip = ir->IntPoint(j);
            eltrans->SetIntPoint(&ip);

            integ_p += ip.weight * eltrans->Weight() * p_gf.GetValue(*eltrans,ip);
            integ_area += ip.weight * eltrans->Weight();
         }
      }
   }          
   
   #ifdef MFEM_USE_MPI
   MPI_Allreduce(MPI_IN_PLACE, &integ_p, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &integ_area, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, vfes->GetComm());
   #endif
   
   MFEM_ASSERT(integ_area > 0.0, "Error: zero area for surface");
   area = integ_area;
   return integ_p / integ_area;
}

//////////////////
//// RESTART /////
//////////////////
void IncompressibleFluidSolver::Restarting(std::string restart_folder_, int step, double &dt, double &time, bool ics)
{
   restart_folder = restart_folder_;
   restarting = ics ? false : true;
   restart_step = step;

   // read in .dat file with time step data and number of processors
   int src_num_procs, step_from_file;
   std::ifstream in_restart;
   std::string restart_dat_file = restart_folder + "/Step_" + std::to_string(step) + "/restart.dat";
   in_restart.open(restart_dat_file.c_str(),std::ios::in);
   if (!in_restart.is_open())
   {
      MFEM_ABORT("Error: unable to open restart file " << restart_dat_file << ". Check the step number is valid.");
   }
   in_restart >> src_num_procs;
   in_restart >> step_from_file;
   in_restart >> time;
   in_restart >> dt;

   // timer info to restore on restart
   int ntimers = 0;
   if (in_restart >> ntimers)
   {
      std::string key;
      double t;
      for (int i = 0; i < ntimers; ++i)
      {
         if (!(in_restart >> key >> t)) { break; }
      }
   }
   in_restart.close();

   MFEM_ASSERT(step_from_file == step, "Error: user provided restart step does not match the step number in the restart file");
}

void IncompressibleFluidSolver::writeRestartFiles(std::string restart_folder_, int restart_save_step_)
{
   restart_folder = restart_folder_;
   restart_save_step = restart_save_step_;
   write_restarts = true;
}

void IncompressibleFluidSolver::writeRestarts(int step, double dt, double time)
{
   sw_restarts.Start();
   // Create directory structure if it doesn't exist
   if (my_rank == 0)
   {      
      std::string step_dir = restart_folder + "/Step_" + std::to_string(step);
   
      try {
         std::filesystem::create_directories(step_dir);
      }
      catch (const std::filesystem::filesystem_error& ex) {
         MFEM_ABORT("Error creating restart directory: " << ex.what());
      }
   }
   #ifdef MFEM_USE_MPI
   MPI_Barrier(vfes->GetComm());
   #endif
   
   // write .dat file with time step data and number of processors
   std::ofstream out_restart;
   std::string restart_dat_file = restart_folder + "/Step_" + std::to_string(step) + "/restart.dat";
   if (my_rank == 0)
   {
      out_restart.open(restart_dat_file.c_str(),std::ios::out);
      if (!out_restart.is_open())
      {
         MFEM_ABORT("Error: unable to open restart file " << restart_dat_file);
      }
   }
   
   int num_procs = 1;
   #ifdef MFEM_USE_MPI
   MPI_Comm_size(vfes->GetComm(), &num_procs);
   #endif

   if (my_rank == 0)
   {
      out_restart << num_procs << std::endl;
      out_restart << step << std::endl;
      out_restart << time << std::endl;
      out_restart << dt << std::endl;

      // Also write current timer values so they can be restored on
      // restart. We (re)register timers here to ensure timer_reg has
      // the up-to-date list and stopwatches.
      RegisterTimers();
      sw_restarts.Stop();
      std::vector<double> timer_times = timer_reg.CollectSeconds();
      sw_restarts.Start();
      const auto &entries = timer_reg.Entries();

      out_restart << entries.size() << std::endl;
      for (size_t i = 0; i < entries.size(); ++i)
      {
         out_restart << entries[i].key << " " << timer_times[i] << std::endl;
      }
      out_restart.close();
   }
   #ifdef MFEM_USE_MPI
   MPI_Barrier(vfes->GetComm());
   #endif

   // Output bc data with history terms (windkessel)
   if (windkessel_bcs.size() > 0)
   {
      std::ofstream out_bc;
      std::string restart_bc_file = restart_folder + "/Step_" + std::to_string(step) + "/windkessel_bcs.dat";
      if (my_rank == 0)
      {
         out_bc.open(restart_bc_file.c_str(), std::ios::out);
         if (!out_bc.is_open())
         {
            MFEM_ABORT("Error: unable to open restart bc data file " << restart_bc_file);
         }

         for (const auto &wbc : windkessel_bcs)
         {
            for (int i = 0; i < wbc.attr.Size(); ++i)
            {
               out_bc << wbc.attr[i] << " ";
            }
            out_bc << std::endl;
            out_bc << wbc.R << " " << wbc.C << " " << wbc.Rd << std::endl;
            out_bc << wbc.hist_terms << std::endl;
            out_bc << wbc.t_nm1 << " " << wbc.t_n << std::endl;
            out_bc << wbc.q_nm1 << " " << wbc.q_n << std::endl;
            out_bc << std::endl; // Blank line as separator
         }
         out_bc.close();
      }
   }

   // Save the mesh and the data
   std::ostringstream oss;
   oss << std::setw(10) << std::setfill('0') << my_rank;
   std::string mname= restart_folder + "/Step_" + std::to_string(step) + "/fluid_"+oss.str()+".mesh";
   std::string vname= restart_folder + "/Step_" + std::to_string(step) + "/velocity_"+oss.str()+".gf";
   std::string pname= restart_folder + "/Step_" + std::to_string(step) + "/pressure_"+oss.str()+".gf";
   std::string vdotname= restart_folder + "/Step_" + std::to_string(step) + "/velocity_dot_"+oss.str()+".gf";
   std::string pdotname= restart_folder + "/Step_" + std::to_string(step) + "/pressure_dot_"+oss.str()+".gf";
   
   std::string solidmname, dname, ddotname, d2dotname;
   if (IFEM)
   {
      solidmname = restart_folder + "/Step_" + std::to_string(step) + "/solid_"+oss.str()+".mesh";
      dname = restart_folder + "/Step_" + std::to_string(step) + "/d_" + oss.str() + ".gf";
      ddotname = restart_folder + "/Step_" + std::to_string(step) + "/ddot_" + oss.str() + ".gf";
      d2dotname = restart_folder + "/Step_" + std::to_string(step) + "/d2dot_" + oss.str() + ".gf";
   }
   std::ofstream sout;

   // Save the mesh
   sout.open(mname.c_str(),std::ios::out);
   sout.precision(20);
   #ifdef MFEM_USE_MPI
   pmesh->ParPrint(sout);
   #else
   pmesh->Print(sout);
   #endif
   sout.close();

   // Save the velocity data
   sout.open(vname.c_str(),std::ios::out);
   sout.precision(20);
   vn_gf.Save(sout);
   sout.close();

   sout.open(vdotname.c_str(),std::ios::out);
   sout.precision(20);
   vdotn_gf.Save(sout);
   sout.close();

   // Save the pressure data
   sout.open(pname.c_str(),std::ios::out);
   sout.precision(20);
   pn_gf.Save(sout);
   sout.close();

   sout.open(pdotname.c_str(),std::ios::out);
   sout.precision(20);
   pdotn_gf.Save(sout);
   sout.close();

   if (IFEM)
   {
      im_solid->writeRestarts(sout, solidmname, dname, ddotname, d2dotname);
   }

   sw_restarts.Start();
}

void IncompressibleFluidSolver::readRestarts(int &step, double &dt, double &time, bool read_bcs)
{
   // this assumes that the fluid mesh is not moving
   // if ALE, then pull back mesh to reference config to map!
   if (ale) MFEM_ABORT("Error: restarting with ALE not implemented yet");

   // Below maps the restart data from potentially different number of processors
   // to the current number of processors   

   // read in .dat file with time step data, number of processors,
   // and (optionally) timer information written by writeRestarts
   int src_num_procs;
   std::ifstream in_restart;
   std::string restart_dat_file = restart_folder + "/Step_" + std::to_string(step) + "/restart.dat";
   in_restart.open(restart_dat_file.c_str(),std::ios::in);
   if (!in_restart.is_open())
   {
      MFEM_ABORT("Error: unable to open restart file " << restart_dat_file);
   }
   in_restart >> src_num_procs;

   // step, time, dt, and then the timer block. Consume these if present so that the stream is
   // positioned correctly before reading timer data.
   int    step_from_file = 0;
   double time_from_file = 0.0;
   double dt_from_file   = 0.0;
   if (in_restart >> step_from_file >> time_from_file >> dt_from_file)
   {
      MFEM_ASSERT(step_from_file == step,
                  "Error: restart step does not match value in restart.dat");
   }

   // Restore per-timer restart times if they were written
   RegisterTimers();
   int ntimers = 0;
   if (in_restart >> ntimers)
   {
      for (int i = 0; i < ntimers; ++i)
      {
         std::string key;
         double t = 0.0;
         if (!(in_restart >> key >> t)) { break; }
         timer_reg.SetRestartTime(key, t);
      }
   }
   in_restart.close();

   std::ifstream in;
   BaseKDTreeNodalProjection* map_v = nullptr;
   BaseKDTreeNodalProjection* map_p = nullptr;
   BaseKDTreeNodalProjection* map_vdotn = nullptr;
   BaseKDTreeNodalProjection* map_pdotn = nullptr;

   BaseKDTreeNodalProjection* map_d = nullptr;
   BaseKDTreeNodalProjection* map_ddot = nullptr;
   BaseKDTreeNodalProjection* map_d2dot = nullptr;

   int dim = pmesh->SpaceDimension();
   if (dim==2)
   {
      map_v = new KDTreeNodalProjection<2>(vn_gf);
      map_p = new KDTreeNodalProjection<2>(pn_gf);
      map_vdotn = new KDTreeNodalProjection<2>(vdotn_gf);
      map_pdotn = new KDTreeNodalProjection<2>(pdotn_gf);

      if (IFEM)
      {
         map_d = new KDTreeNodalProjection<2>(*im_solid->GetCurrentDisplacement());
         map_ddot = new KDTreeNodalProjection<2>(*im_solid->GetCurrentVelocity());
         map_d2dot = new KDTreeNodalProjection<2>(*im_solid->GetCurrentAcceleration());
      }
   }
   else
   {
      map_v = new KDTreeNodalProjection<3>(vn_gf);
      map_p = new KDTreeNodalProjection<3>(pn_gf);
      map_vdotn = new KDTreeNodalProjection<3>(vdotn_gf);
      map_pdotn = new KDTreeNodalProjection<3>(pdotn_gf);

      if (IFEM)
      {
         map_d = new KDTreeNodalProjection<3>(*im_solid->GetCurrentDisplacement());
         map_ddot = new KDTreeNodalProjection<3>(*im_solid->GetCurrentVelocity());
         map_d2dot = new KDTreeNodalProjection<3>(*im_solid->GetCurrentAcceleration());
      }
   }

   // loop over all source processors and read in data, project to current mesh
   for (int p=0; p<src_num_procs; p++)
   {
      std::ostringstream oss;
      oss << std::setw(10) << std::setfill('0') << p;
      std::string mname= restart_folder + "/Step_" + std::to_string(step) + "/fluid_"+oss.str()+".mesh";
      std::string vname= restart_folder + "/Step_" + std::to_string(step) + "/velocity_"+oss.str()+".gf";
      std::string pname= restart_folder + "/Step_" + std::to_string(step) + "/pressure_"+oss.str()+".gf";
      std::string vdotname= restart_folder + "/Step_" + std::to_string(step) + "/velocity_dot_"+oss.str()+".gf";
      std::string pdotname= restart_folder + "/Step_" + std::to_string(step) + "/pressure_dot_"+oss.str()+".gf";

      std::string solidmname, dname, ddotname, d2dotname;
      if (IFEM)
      {
         solidmname = restart_folder + "/Step_" + std::to_string(step) + "/solid_"+oss.str()+".mesh";
         dname = restart_folder + "/Step_" + std::to_string(step) + "/d_" + oss.str() + ".gf";
         ddotname = restart_folder + "/Step_" + std::to_string(step) + "/ddot_" + oss.str() + ".gf";
         d2dotname = restart_folder + "/Step_" + std::to_string(step) + "/d2dot_" + oss.str() + ".gf";
      }

      // Read the mesh
      Mesh lmesh;
      in.open(mname.c_str(),std::ios::in);
      lmesh.Load(in);
      in.close();

      in.open(vname.c_str(),std::ios::in);
      GridFunction v_in_gf(&lmesh,in);
      in.close();

      in.open(vdotname.c_str(),std::ios::in);
      GridFunction vdot_in_gf(&lmesh,in);
      in.close();

      in.open(pname.c_str(),std::ios::in);
      GridFunction p_in_gf(&lmesh,in);
      in.close();

      in.open(pdotname.c_str(),std::ios::in);
      GridFunction pdot_in_gf(&lmesh,in);
      in.close();

      // Project the grid functions
      double tol = 1e-10;
      map_v->Project(v_in_gf,tol);
      map_p->Project(p_in_gf,tol);
      map_vdotn->Project(vdot_in_gf,tol);
      map_pdotn->Project(pdot_in_gf,tol);

      if (IFEM)
      {
         im_solid->readRestarts(in, map_d, map_ddot, map_d2dot, solidmname, dname, ddotname, d2dotname);
      }
   }
   delete map_v;
   delete map_p;
   delete map_vdotn;
   delete map_pdotn;
   if (IFEM)
   {
      delete map_d;
      delete map_ddot;
      delete map_d2dot;
   }

   #ifdef MFEM_USE_MPI
   MPI_Barrier(vfes->GetComm());
   #endif

   // update yn and ydotn
   Vector &vn = yn.GetBlock(0);
   Vector &pn = yn.GetBlock(1);
   #ifdef MFEM_USE_MPI
   vn_gf.ParallelProject(vn);
   pn_gf.ParallelProject(pn);
   #else
   vn_gf.GetTrueDofs(vn);
   pn_gf.GetTrueDofs(pn);
   #endif

   Vector &vdotn = ydotn.GetBlock(0);
   Vector &pdotn = ydotn.GetBlock(1);
   #ifdef MFEM_USE_MPI
   vdotn_gf.ParallelProject(vdotn);
   pdotn_gf.ParallelProject(pdotn);
   #else
   vdotn_gf.GetTrueDofs(vdotn);
   pdotn_gf.GetTrueDofs(pdotn);
   #endif

   if (IFEM)
   {
      im_solid->UpdateSolidStateFromRestart();
   }

   // read bc data
   if (windkessel_bcs.size() > 0 && read_bcs)
   {
      std::string restart_bc_file = restart_folder + "/Step_" + std::to_string(step) + "/windkessel_bcs.dat";
      std::ifstream in_bc(restart_bc_file.c_str());
      if (!in_bc.is_open())
      {
         MFEM_ABORT("Error: unable to open restart bc data file " << restart_bc_file);
      }

      std::string line;
      while (std::getline(in_bc, line))
      {
         if (line.empty()) continue; // Skip blank lines

         // Read attr array (all ints on one line)
         Array<int> attr;
         std::istringstream iss(line);
         int val;
         while (iss >> val)
         {
            attr.Append(val);
         }

         // Next line: R C Rd
         std::getline(in_bc, line);
         double R = -1.0, C = -1.0, Rd = -1.0;
         std::istringstream iss_params(line);
         iss_params >> R >> C >> Rd;
         if (R < 0 || C < 0 || Rd < 0)
         {
            MFEM_ABORT("Error reading windkessel parameters from restart bc file");
         }

         std::getline(in_bc, line); // hist_terms
         double hist_terms = std::stod(line);

         std::getline(in_bc, line); // times
         double t_nm1, t_n;
         std::istringstream iss_time(line);
         iss_time >> t_nm1 >> t_n;

         std::getline(in_bc, line); // flow rates
         double q_nm1, q_n;
         std::istringstream iss_flow(line);
         iss_flow >> q_nm1 >> q_n;

         // Update the matching windkessel bc
         // Note: assumes unique (attr, Rp, C, Rd) tuple
         for (auto &wbc : windkessel_bcs)
         {
            // check if attributes and parameters match
            bool attr_match = true;
            if (attr.Size() != wbc.attr.Size()) attr_match = false;
            for (int i = 0; i < attr.Size(); i++)
               if (attr[i] != wbc.attr[i]) attr_match = false;
            if (!attr_match) continue;
            if (wbc.R != R || wbc.C != C || wbc.Rd != Rd) continue;

            // update
            wbc.hist_terms = hist_terms;
            wbc.t_nm1 = t_nm1; wbc.t_n = t_n;
            wbc.q_nm1 = q_nm1; wbc.q_n = q_n;

            if (ale)
            {
                  wbc.AssemblePVEC();
                  wbc.AssembleSVEC();
            }
            break;
         }
         std::getline(in_bc, line); // blank line separator
      }
      in_bc.close();
   }
}

//////////////////
/// DESTRUCTOR ///
//////////////////

IncompressibleFluidSolver::~IncompressibleFluidSolver()
{
   delete vfec;
   delete pfec;
   delete vfes;
   delete pfes;
   delete mass_lf;
   
   for (int i=0; i<fes_stab.Size(); i++)
   {
      fes_stab[i] = nullptr;
      if (!total_integ) total_ess_tdofs[i] = nullptr;
   }

   delete f_form;
   delete L_L_form;
   
   delete res;
   delete newton_solver;
   delete solver;
   delete solver_monitor;
   delete solver_wrap;
   delete jacobian_prec;

   delete cvbc_op;
   delete cvbc_PC;

   delete FluidNF;

   delete im_solid;
}