#include "predictor_corrector_newton_fluid.hpp"

using namespace mfem;
using namespace std;
using namespace ins_solver_serial;
using namespace immersed_solid_serial;

void PC_NewtonSolverFluid::Mult(const Vector &b, Vector &x) const
{
   MFEM_ABORT("not implemented");
}

void PC_NewtonSolverFluid::NewtonSolve(const double &dt, const double &time) const
{
   MFEM_ASSERT(oper != NULL, "the Operator is not set (use SetOperator).");
   MFEM_ASSERT(prec != NULL, "the Solver is not set (use SetSolver).");

   prec->iterative_mode = false;

   // Fix Cf operator at t_n
   if (IFEM)
   {
      im_solid->ComputeMortarIntegral();
   }

   // x_{i+1} = x_i - [DF(x_i)]^{-1} [F(x_i)-b]
   for (it = 0; true; it++)
   {
      // initiator step
      Initiator();

      // Set weak bcs
      // if mesh is not moving, only set on it = 0
      // to t_n+alpha_f = t_n + alpha_f*dt, time = t_n+1
      // t_n+alpha_f = time - dt + alpha_f*dt
      if (it == 0) SetWeakBCs(time - dt + alpha_f*dt,dt);    

      // Form Linear System
      // for gen alpha at midpoint
      sw_res->Start();
      FormRHS(*yg, *ydotg, dt, r); // minus sign NOT included in this term!   
      sw_res->Stop();

      // check residual before solving
      bool temp = CheckConvergence();
      PrintNewtonSolveInfo();
      if (temp)
      {
         if (print_newton)
         {
            #ifdef MFEM_USE_MPI
            if (my_rank == 0) mfem::out << "\n";
            #else
            mfem::out << "\n";
            #endif
         }
         break;
      }       
      
      // form LHS for gen alpha at midpoint
      sw_tan->Start();
      FormLHS(*yg, *ydotg, dt);
      sw_tan->Stop();

      // Solve
      sw_linsolve->Start();
      SolveLinearSystem();
      sw_linsolve->Stop();
      PrintLinSolveInfo();

      // Corrector
      Correctors(dt,time);
   }

   final_iter = it;
   final_norm = norm;

   if (!converged)
   {
      if (no_error_nonconvergence)
      {
         #ifdef MFEM_USE_MPI
         if (my_rank == 0) mfem::out << "Newton: No convergence!\n";
         #else
         mfem::out << "Newton: No convergence!\n";
         #endif
      }
      else
      {
         #ifdef MFEM_USE_MPI
         if (my_rank == 0) MFEM_ABORT("Newton: No convergence!\n");
         #else
         MFEM_ABORT("Newton: No convergence!\n");
         #endif
      }
   }
}

void PC_NewtonSolverFluid::PrintNorms() const
{  
   rblock.Update(ynp1->GetData(),block_trueOffsets);
   Vector &unp1 = rblock.GetBlock(0);
   Vector &pnp1 = rblock.GetBlock(1);
   #ifdef MFEM_USE_MPI
   double glob_vel_norm = mfem::InnerProduct(MPI_COMM_WORLD,unp1,unp1);
   double glob_p_norm = mfem::InnerProduct(MPI_COMM_WORLD,pnp1,pnp1);

   if (Mpi::Root())
   {
      std::cout << setprecision(16);
      // the first is not the L2 norm!!!
      std::cout << "|| u_h || = " << glob_vel_norm << "\n";
      std::cout << "|| p_h || = " << glob_p_norm << "\n";
   }
   #else
   double glob_vel_norm = mfem::InnerProduct(unp1,unp1);
   double glob_p_norm = mfem::InnerProduct(pnp1,pnp1);

   std::cout << setprecision(16);
   // the first is not the L2 norm!!!
   std::cout << "|| u_h || = " << glob_vel_norm << "\n";
   std::cout << "|| p_h || = " << glob_p_norm << "\n";
   #endif   
}

void PC_NewtonSolverFluid::Initiator() const
{
   // compute midpoint values
   add((1.0 - alpha_f),*yn,alpha_f,*ynp1,*yg);
   add((1.0 - alpha_m),*ydotn,alpha_m,*ydotnp1,*ydotg);

   if(IFEM)
   {
      im_solid->Initiator();
   }
}

void PC_NewtonSolverFluid::SetWeakBCs(const double &time, double dt) const
{
   // traction bcs
   for (auto &vel_nbc : *vel_nbcs)
   {
      if (vel_nbc.coeff)
      {
         vel_nbc.coeff->SetTime(time);
      }
      else if(vel_nbc.mcoeff)
      {
         vel_nbc.mcoeff->SetTime(time);
      }
      else if (vel_nbc.scalar_coeff)
      {
         vel_nbc.scalar_coeff->SetTime(time);
      }
      else
      {
         MFEM_ABORT("Error in defining velocity neumann bcs");
      }      
   }
   L_L_form->Assemble();
   #ifdef MFEM_USE_MPI
   dynamic_cast<ParLinearForm *>(L_L_form)->ParallelAssemble(*ft);
   #else
   ft->Set(1.0,*L_L_form);
   #endif

   // body forces
   for (auto &accel_term : *accel_terms)
   {
      accel_term.coeff->SetTime(time);
   }
   if (!total_integ)
   {      
      f_form->Assemble();
      #ifdef MFEM_USE_MPI
      dynamic_cast<ParLinearForm *>(f_form)->ParallelAssemble(*fb);
      #else
      fb->Set(1.0,*f_form);
      #endif
   }   

   // Cardiovascular bcs
   if (resistance_bc_flag)
   {
      // nothing needed here
   }
   if (windkessel_bc_flag)
   {
      for (auto &wbc : *windkessel_bcs)
      {
         // set times
         wbc.dt = dt;

         // time = t_n + alpha_f*dt
         wbc.t_npalphaf = time;
         wbc.t_np1 = time - alpha_f*dt + dt; // t_n+1 = t_n + dt

         // set integ hist terms (multiply by e^(-t_n+alpha_F / (Rd*C)))
         wbc.integ_hist_terms = wbc.hist_terms*exp(-wbc.t_npalphaf / (wbc.Rd*wbc.C) );
      }
   }

   // immersed solid
   if (IFEM)
   {
      im_solid->SetWeakBCs(time);
   }
}

void PC_NewtonSolverFluid::FormRHS(Vector &y, Vector &ydot, double dt, Vector &res) const
{
   double am = 1.0;
   res_op->Mult_PC(y,ydot,*ft,*fb,am,res);
}

void PC_NewtonSolverFluid::FormLHS(Vector &y, Vector &ydot, double dt) const
{
   // Tangent Stiffness //
   double am = alpha_m;
   double af = alpha_f*gamma*dt;

   grad = &res_op->GetGradient_PC(y,ydot,dt,am,af);
   // reset linear solver operator (which is prec here)
   prec->SetOperator(*grad);

   // // setup pre-preconditioning 
   // if (jacobi_pre_precond)
   // {
   //    #ifdefdef MFEM_USE_MPI
   //    MFEM_ASSERT(direct_solver == false, "Jacobi pre-preconditioner only works with iterative solvers");
   //    BlockOperator *block_grad = dynamic_cast<BlockOperator *>(grad);
   //    if (block_grad)
   //    {
   //       HypreParMatrix *A_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(0,0)));
   //       HypreParMatrix *G_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(0,1)));
   //       HypreParMatrix *D_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(1,0)));
   //       HypreParMatrix *C_stab = nullptr;
   //       if (equal_order) C_stab = const_cast<HypreParMatrix *>(dynamic_cast<const HypreParMatrix *>(&grad_stab->GetBlock(1,1)));
         
   //       // Get diagonal of A and C
   //       Vector &Ad = res_op->GetADiagonal();
   //       Vector &Cd = res_op->GetCDiagonal();

   //       // Scale matrices
   //       HypreParMatrixInvScaleRows(A_stab,Ad);
   //       if (equal_order) HypreParMatrixInvScaleRows(C_stab,Cd);

   //       HypreParMatrix Wm()

   //    }
   //    else
   //    {
   //       MFEM_ABORT("Jacobi pre-preconditioner only works with BlockOperator");
   //    }
   //    #else
   //    MFEM_ABORT("TODO");
   //    #endif
   // }
}

void PC_NewtonSolverFluid::SolveLinearSystem() const
{
   prec->Mult(r,c);
   if (!direct_solver)
   {
      if (dynamic_cast<IterativeSolver *>(prec))
      {
         final_it_ls_fs = dynamic_cast<IterativeSolver *>(prec)->GetNumIterations();
         final_atol_ls_fs = dynamic_cast<IterativeSolver *>(prec)->GetFinalNorm();
         final_rtol_ls_fs = dynamic_cast<IterativeSolver *>(prec)->GetFinalRelNorm(); 
      }
      else
      {
         final_it_ls_fs = dynamic_cast<BlockOrthoSolver*>(prec)->GetNumIterations();
         final_atol_ls_fs = dynamic_cast<BlockOrthoSolver*>(prec)->GetFinalNorm();
         final_rtol_ls_fs = dynamic_cast<BlockOrthoSolver*>(prec)->GetFinalRelNorm(); 
      }
   } 
   else
   {
      final_it_ls_fs = 0;
      final_atol_ls_fs = 0;
      final_rtol_ls_fs = 0;
   }
}

void PC_NewtonSolverFluid::Correctors(const double &dt, const double &time) const
{
   if (IFEM)
   {
      delblock.Update(c.GetData(),block_trueOffsets);
      Vector &delv = delblock.GetBlock(0);

      // project update to solid mesh and update quantities
      im_solid->Correctors(delv);
   }

   ydotnp1->Add(-1.0,c);
   ynp1->Add(-gamma*dt,c);
}

bool PC_NewtonSolverFluid::CheckConvergence() const
{
   norm = Norm(r);

   if (it == 0)
   {
      norm0 = initial_norm = norm;
      norm_goal = std::max(rel_tol*norm, abs_tol);
   } 
   MFEM_ASSERT(IsFinite(norm), "norm = " << norm);
   if (print_options.iterations)
   {
      mfem::out << "Newton iteration " << setw(2) << it
                  << " : ||r|| = " << norm;
      if (it > 0)
      {
         mfem::out << ", ||r||/||r_0|| = " << norm/norm0;
      }
      mfem::out << '\n';
   }
   if (norm <= norm_goal && it >= min_iter)
   {
      converged = true;
      return true;
      // break;
   }
   if (it >= max_iter)
   {
      converged = false;
      return true;
      // break;
   }
   return false;
}

void PC_NewtonSolverFluid::PrintNewtonSolveInfo() const
{
   #ifdef MFEM_USE_MPI
   if (print_newton && my_rank == 0)
   #else
   if (print_newton)
   #endif
   {   
      if (it == 0)
      {
         mfem::out << std::endl;
         mfem::out << std::setw(7) << "Eqn  " << std::setw(3) << "It" << std::setw(9)
                << "AbsRes" << std::setw(11) << "RelRes" << std::setw(9) << "LS_It" 
                << std::setw(11) << "LS_AbsRes" << std::setw(11) << "LS_RelRes"
                << "\n";
      }

      rel_norm = std::isnan(norm/norm0) ? -1.0 : norm/norm0;
      if(std::isinf(rel_norm)) rel_norm = -1.0;

      mfem::out << std::setw(5) << "Fluid" << std::setw(5) << std::fixed
                << it << "   " << std::setw(3) << std::setprecision(2)
                << std::scientific << norm << "   " << rel_norm; //<< std::flush;

      mfem::out << std::setprecision(8);
   }
}

void PC_NewtonSolverFluid::PrintLinSolveInfo() const
{
   #ifdef MFEM_USE_MPI
   if (print_newton && my_rank == 0)
   #else
   if (print_newton)
   #endif
   {
      rel_norm_ls = std::isnan(final_rtol_ls_fs) ? -1.0 : final_rtol_ls_fs;
      mfem::out << "  " << std::setw(5) << std::fixed
                     << final_it_ls_fs << "  " << std::setw(3) << std::setprecision(2)
                     << std::scientific << final_atol_ls_fs << "   " << rel_norm_ls 
                     << "\n"; //<< std::flush;
      mfem::out << std::setprecision(8);
   }
}