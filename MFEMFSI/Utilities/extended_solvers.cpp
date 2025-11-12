#include "extended_solvers.hpp"
#include <general/forall.hpp>
#include <linalg/solvers.hpp>
#include <iomanip>

using namespace mfem;
using namespace std;
using namespace modified_gmres_serial;

void ExtendedNewtonSolver::Mult(const Vector &b, Vector &x) const
{
   MFEM_ASSERT(oper != NULL, "the Operator is not set (use SetOperator).");
   MFEM_ASSERT(prec != NULL, "the Solver is not set (use SetSolver).");

   int it;
   double norm0, norm, norm_goal;
   const bool have_b = (b.Size() == Height());

   if (!iterative_mode)
   {
      x = 0.0;
   }

   ProcessNewState(x);

   oper->Mult(x, r);
   if (have_b)
   {
      r -= b;
   }

   // std::cout << setprecision(16);
   // std::cout << "Initial residual :" << std::endl;
   // r.Print(cout);
   // std::cout << "---------------------\n";

   norm0 = norm = initial_norm = Norm(r);
   if (print_options.first_and_last && !print_options.iterations)
   {
      mfem::out << "Newton iteration " << setw(2) << 0
                << " : ||r|| = " << norm << "...\n";
   }
   norm_goal = std::max(rel_tol*norm, abs_tol);

   prec->iterative_mode = false;

   // x_{i+1} = x_i - [DF(x_i)]^{-1} [F(x_i)-b]
   for (it = 0; true; it++)
   {
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
      Monitor(it, norm, r, x);

      if (apply_min_iter_after1)
      {
         if (norm <= norm_goal)
         {
            converged = true;
            break;
         }
         apply_min_iter_after1 = false;
      }
      else
      {
         if (norm <= norm_goal && it >= min_iter)
         {
            converged = true;
            break;
         }
      }      

      if (it >= max_iter)
      {
         converged = false;
         break;
      }

      // MFEM_PERF_BEGIN("Form Jacobian");
      grad = &oper->GetGradient(x);
      // MFEM_PERF_END("Form Jacobian");

      // // printing out the matrix
      // if (it > 0)
      // {
      //    string file_prefix = "my_mat";
      //    dynamic_cast<HypreParMatrix *>(grad)->Print(file_prefix.c_str());
      // }
      // string file_prefix = "my_mat";
      // dynamic_cast<HypreParMatrix *>(grad)->Print(file_prefix.c_str());

      // grad->PrintMatlab(cout);
      

      prec->SetOperator(*grad);

      if (lin_rtol_type)
      {
         AdaptiveLinRtolPreSolve(x, it, norm);
      }

      // MFEM_PERF_BEGIN("Solve Linearized Equation");
      prec->Mult(r, c); // c = [DF(x_i)]^{-1} [F(x_i)-b]
      // MFEM_PERF_END("Solve Linearized Equation");

      if (lin_rtol_type)
      {
         AdaptiveLinRtolPostSolve(c, r, it, norm);
      }

      const double c_scale = ComputeScalingFactor(x, b);
      if (c_scale == 0.0)
      {
         converged = false;
         break;
      }

      // mfem::out << setprecision(16);
      // mfem::out << "before update:" << std::endl;
      // x.Print();

      // mfem::out << "delta:" << std::endl;
      // c.Print();




      add(x, -c_scale, c, x);


      // mfem::out << "after update:" << std::endl;
      // x.Print();


      ProcessNewState(x);
      if (gen_alpha)
      {
         UpdateGenAlphaVelocity(c); // pass in delta [udot, p]
      }

      // MFEM_PERF_BEGIN("Compute RHS");
      oper->Mult(x, r);
      // MFEM_PERF_END("Compute RHS");

      // mfem::out << setprecision(16);
      // mfem::out << "after update residual:" << std::endl;
      // r.Print(mfem::out);

      if (have_b)
      {
         r -= b;
      }
      norm = Norm(r);
   }

   final_iter = it;
   final_norm = norm;

   if (print_options.summary || (!converged && print_options.warnings) ||
       print_options.first_and_last)
   {
      mfem::out << "Newton: Number of iterations: " << final_iter << '\n'
                << "   ||r|| = " << final_norm << '\n';
   }
   if (!converged && (print_options.summary || print_options.warnings))
   {
      mfem::out << "Newton: No convergence!\n";
   }
}

inline void GeneratePlaneRotation(double &dx, double &dy,
                                  double &cs, double &sn)
{
   if (dy == 0.0)
   {
      cs = 1.0;
      sn = 0.0;
   }
   else if (fabs(dy) > fabs(dx))
   {
      double temp = dx / dy;
      sn = 1.0 / sqrt( 1.0 + temp*temp );
      cs = temp * sn;
   }
   else
   {
      double temp = dy / dx;
      cs = 1.0 / sqrt( 1.0 + temp*temp );
      sn = temp * cs;
   }
}

inline void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
   double temp = cs * dx + sn * dy;
   dy = -sn * dx + cs * dy;
   dx = temp;
}

inline void Update(Vector &x, int k, DenseMatrix &h, Vector &s,
                   Array<Vector*> &v)
{
   Vector y(s);

   // Backsolve:
   for (int i = k; i >= 0; i--)
   {
      y(i) /= h(i,i);
      for (int j = i - 1; j >= 0; j--)
      {
         y(j) -= h(j,i) * y(i);
      }
   }

   for (int j = 0; j <= k; j++)
   {
      x.Add(y(j), *v[j]);
   }
}

void ExtendedFGMRESSolver::Mult(const Vector &b, Vector &x) const
{
   DenseMatrix H(m+1,m);
   Vector s(m+1), cs(m+1), sn(m+1);
   Vector r(b.Size());

   int i, j, k;


   if (iterative_mode)
   {
      oper->Mult(x, r);
      subtract(b,r,r);
   }
   else
   {
      x = 0.;
      r = b;
   }
   double beta = initial_norm = Norm(r);  // beta = ||r||
   // We need to preallocate this to report the correct result in the case of
   // no convergence.
   double resid;
   MFEM_ASSERT(IsFinite(beta), "beta = " << beta);

   final_norm = std::max(rel_tol*beta, abs_tol);

   // if (beta <= final_norm)
   // {
   //    final_norm = beta;
   //    final_iter = 0;
   //    converged = true;
   //    return;
   // }

   if (print_options.iterations || print_options.first_and_last)
   {
      mfem::out << "   Pass : " << setw(2) << 1
                << "   Iteration : " << setw(3) << 0
                << "  || r || = " << beta
                << (print_options.first_and_last ? " ...\n" : "\n");
   }

   Monitor(0, beta, r, x);

   Array<Vector*> v(m+1);
   Array<Vector*> z(m+1);
   for (i= 0; i<=m; i++)
   {
      v[i] = NULL;
      z[i] = NULL;
   }

   j = 1;
   while (j <= max_iter)
   {
      if (v[0] == NULL) { v[0] = new Vector(b.Size()); }
      (*v[0]) = 0.0;
      v[0] -> Add (1.0/beta, r);   // v[0] = r / ||r||
      s = 0.0; s(0) = beta;

      for (i = 0; i < m && j <= max_iter; i++, j++)
      {

         if (z[i] == NULL) { z[i] = new Vector(b.Size()); }
         (*z[i]) = 0.0;

         if (prec)
         {
            // MFEM_PERF_BEGIN("Preconditioner Apply");
            prec->Mult(*v[i], *z[i]);
            // MFEM_PERF_END("Preconditioner Apply");
         }
         else
         {
            (*z[i]) = (*v[i]);
         }
         // MFEM_PERF_BEGIN("Jacobian Apply");
         oper->Mult(*z[i], r);
         // MFEM_PERF_END("Jacobian Apply");

         for (k = 0; k <= i; k++)
         {
            H(k,i) = Dot( r, *v[k]); // H(k,i) = r * v[k]
            r.Add(-H(k,i), (*v[k])); // r -= H(k,i) * v[k]
         }

         H(i+1,i)  = Norm(r);       // H(i+1,i) = ||r||
         if (v[i+1] == NULL) { v[i+1] = new Vector(b.Size()); }
         (*v[i+1]) = 0.0;
         v[i+1] -> Add (1.0/H(i+1,i), r); // v[i+1] = r / H(i+1,i)

         for (k = 0; k < i; k++)
         {
            ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
         }

         GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
         ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
         ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

         resid = fabs(s(i+1));
         MFEM_ASSERT(IsFinite(resid), "resid = " << resid);
         if (print_options.iterations || (print_options.first_and_last &&
                                          resid <= final_norm))
         {
            mfem::out << "   Pass : " << setw(2) << (j-1)/m+1
                      << "   Iteration : " << setw(3) << j
                      << "  || r || = " << resid << endl;
         }
         Monitor(j, resid, r, x, resid <= final_norm);

         if (resid <= final_norm && j >= min_iter)
         {
            Update(x, i, H, s, z);
            final_norm = resid;
            final_iter = j;
            converged = true;

            if (print_options.summary)
            {
               mfem::out << "FGMRES: Number of iterations: " << final_iter << '\n';
            }

            for (i= 0; i<=m; i++)
            {
               if (v[i]) { delete v[i]; }
               if (z[i]) { delete z[i]; }
            }
            return;
         }
      }

      if (print_options.iterations)
      {
         mfem::out << "Restarting..." << endl;
      }

      Update(x, i-1, H, s, z);

      oper->Mult(x, r);
      subtract(b,r,r);
      beta = Norm(r);
      MFEM_ASSERT(IsFinite(beta), "beta = " << beta);
      if (beta <= final_norm && j >= min_iter)
      {
         final_norm = beta;
         final_iter = j;
         converged = true;

         break;
      }
   }

   if (converged == false)
   {
      final_norm = beta;
      final_iter = -1;
   }

   // Clean buffers up
   for (i = 0; i <= m; i++)
   {
      if (v[i]) { delete v[i]; }
      if (z[i]) { delete z[i]; }
   }
   converged = false;

   // Note: j is off by one when we arrive here
   if (!print_options.iterations && print_options.first_and_last)
   {
      mfem::out << "   Pass : " << setw(2) << (j-1)/m+1
                << "   Iteration : " << setw(3) << j-1
                << "  || r || = " << resid << endl;
   }
   if (print_options.summary || (print_options.warnings && !converged))
   {
      mfem::out << "FGMRES: Number of iterations: " << j-1 << '\n';
   }
   if (print_options.warnings && !converged)
   {
      mfem::out << "FGMRES: No convergence!\n";
   }
}

void ExtendedGMRESSolver::Mult(const Vector &b, Vector &x) const
{
   // Generalized Minimum Residual method following the algorithm
   // on p. 20 of the SIAM Templates book.

   int n = width;

   DenseMatrix H(m+1, m);
   Vector s(m+1), cs(m+1), sn(m+1);
   Vector r(n), w(n);
   Array<Vector *> v;

   double resid;
   int i, j, k;

   if (iterative_mode)
   {
      oper->Mult(x, r);
   }
   else
   {
      x = 0.0;
   }

   if (prec)
   {
      if (iterative_mode)
      {
         subtract(b, r, w);
         prec->Mult(w, r);    // r = M (b - A x)
      }
      else
      {
         prec->Mult(b, r);
      }
   }
   else
   {
      if (iterative_mode)
      {
         subtract(b, r, r);
      }
      else
      {
         r = b;
      }
   }
   double beta = initial_norm = Norm(r);  // beta = ||r||
   MFEM_ASSERT(IsFinite(beta), "beta = " << beta);

   final_norm = std::max(rel_tol*beta, abs_tol);

   // if (beta <= final_norm)
   // {
   //    final_norm = beta;
   //    final_iter = 0;
   //    converged = true;
   //    j = 0;
   //    resid = beta;
   //    goto finish;
   // }

   if (print_options.iterations || print_options.first_and_last)
   {
      mfem::out << "   Pass : " << setw(2) << 1
                << "   Iteration : " << setw(3) << 0
                << "  ||B r|| = " << beta
                << (print_options.first_and_last ? " ...\n" : "\n");
   }

   Monitor(0, beta, r, x);

   v.SetSize(m+1, NULL);

   for (j = 1; j <= max_iter; )
   {
      if (v[0] == NULL) { v[0] = new Vector(n); }
      v[0]->Set(1.0/beta, r);
      s = 0.0; s(0) = beta;

      for (i = 0; i < m && j <= max_iter; i++, j++)
      {
         if (prec)
         {
            // MFEM_PERF_BEGIN("Jacobian Apply 1");
            oper->Mult(*v[i], r);
            // MFEM_PERF_END("Jacobian Apply 1");
            // MFEM_PERF_BEGIN("Preconditioner Apply 1");
            prec->Mult(r, w);        // w = M A v[i]
            // MFEM_PERF_END("Preconditioner Apply 1");
         }
         else
         {
            oper->Mult(*v[i], w);
         }

         for (k = 0; k <= i; k++)
         {
            H(k,i) = Dot(w, *v[k]);  // H(k,i) = w * v[k]
            w.Add(-H(k,i), *v[k]);   // w -= H(k,i) * v[k]
         }

         H(i+1,i) = Norm(w);           // H(i+1,i) = ||w||
         MFEM_ASSERT(IsFinite(H(i+1,i)), "Norm(w) = " << H(i+1,i));
         if (v[i+1] == NULL) { v[i+1] = new Vector(n); }
         v[i+1]->Set(1.0/H(i+1,i), w); // v[i+1] = w / H(i+1,i)

         for (k = 0; k < i; k++)
         {
            ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
         }

         GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
         ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
         ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

         resid = fabs(s(i+1));
         MFEM_ASSERT(IsFinite(resid), "resid = " << resid);

         if (resid <= final_norm && j >= min_iter)
         {
            Update(x, i, H, s, v);
            final_norm = resid;
            final_iter = j;
            converged = true;
            goto finish;
         }

         if (print_options.iterations)
         {
            mfem::out << "   Pass : " << setw(2) << (j-1)/m+1
                      << "   Iteration : " << setw(3) << j
                      << "  ||B r|| = " << resid << '\n';
         }

         Monitor(j, resid, r, x);
      }

      if (print_options.iterations && j <= max_iter)
      {
         mfem::out << "Restarting..." << '\n';
      }

      Update(x, i-1, H, s, v);

      // MFEM_PERF_BEGIN("Jacobian Apply 2");
      oper->Mult(x, r);
      // MFEM_PERF_END("Jacobian Apply 2");
      if (prec)
      {
         subtract(b, r, w);
         // MFEM_PERF_BEGIN("Preconditioner Apply 2");
         prec->Mult(w, r);    // r = M (b - A x)
         // MFEM_PERF_END("Preconditioner Apply 2");
      }
      else
      {
         subtract(b, r, r);
      }
      beta = Norm(r);         // beta = ||r||
      MFEM_ASSERT(IsFinite(beta), "beta = " << beta);
      if (beta <= final_norm && j >= min_iter)
      {
         final_norm = beta;
         final_iter = j;
         converged = true;
         goto finish;
      }
   }

   final_norm = beta;
   final_iter = max_iter;
   converged = false;

finish:
   if ((print_options.iterations && converged) || print_options.first_and_last)
   {
      mfem::out << "   Pass : " << setw(2) << (j-1)/m+1
                << "   Iteration : " << setw(3) << final_iter
                << "  ||B r|| = " << resid << '\n';
   }
   if (print_options.summary || (print_options.warnings && !converged))
   {
      mfem::out << "GMRES: Number of iterations: " << final_iter << '\n';
   }
   if (print_options.warnings && !converged)
   {
      mfem::out << "GMRES: No convergence!\n";
   }

   Monitor(final_iter, final_norm, r, x, true);

   for (i = 0; i < v.Size(); i++)
   {
      delete v[i];
   }
}

void ModifiedGMRESSolver::Mult(const Vector &b, Vector &x) const
{
   // Generalized Minimum Residual method following the algorithm
   // on p. 20 of the SIAM Templates book.

   int n = width;

   DenseMatrix H(m+1, m);
   Vector s(m+1), cs(m+1), sn(m+1);
   Vector r(n), w(n);
   Array<Vector *> v;

   int i, j, k;

   if (iterative_mode)
   {
      oper->Mult(x, r);
   }
   else
   {
      x = 0.0;
   }

   if (prec)
   {
      if (iterative_mode)
      {
         subtract(b, r, w);
         prec->Mult(w, r);    // r = M (b - A x)
      }
      else
      {
         prec->Mult(b, r);
      }
   }
   else
   {
      if (iterative_mode)
      {
         subtract(b, r, r);
      }
      else
      {
         r = b;
      }
   }
   double beta = initial_norm = Norm(r);  // beta = ||r||
   MFEM_ASSERT(IsFinite(beta), "beta = " << beta);

   final_norm = std::max(rel_tol*beta, abs_tol);

   if (beta <= final_norm)
   {
      final_norm = beta;
      final_iter = 0;
      converged = true;
      j = 0;
      goto finish;
   }

   if (print_options.iterations || print_options.first_and_last)
   {
      mfem::out << "   Pass : " << setw(2) << 1
               << "   Iteration : " << setw(3) << 0
               << "  ||B r|| = " << beta
               << (print_options.first_and_last ? " ...\n" : "\n");
   }

   Monitor(0, beta, r, x);

   v.SetSize(m+1, NULL);

   for (j = 1; j <= max_iter; )
   {
      if (v[0] == NULL) { v[0] = new Vector(n); }
      v[0]->Set(1.0/beta, r);
      s = 0.0; s(0) = beta;

      for (i = 0; i < m && j <= max_iter; i++, j++)
      {
         if (prec)
         {
            oper->Mult(*v[i], r);
            prec->Mult(r, w);        // w = M A v[i]
         }
         else
         {
            oper->Mult(*v[i], w);
         }

         for (k = 0; k <= i; k++)
         {
            H(k,i) = Dot(w, *v[k]);  // H(k,i) = w * v[k]
            w.Add(-H(k,i), *v[k]);   // w -= H(k,i) * v[k]
         }

         H(i+1,i) = Norm(w);           // H(i+1,i) = ||w||
         MFEM_ASSERT(IsFinite(H(i+1,i)), "Norm(w) = " << H(i+1,i));
         if (v[i+1] == NULL) { v[i+1] = new Vector(n); }
         v[i+1]->Set(1.0/H(i+1,i), w); // v[i+1] = w / H(i+1,i)

         for (k = 0; k < i; k++)
         {
            ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
         }

         GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
         ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
         ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

         const double resid = fabs(s(i+1));
         MFEM_ASSERT(IsFinite(resid), "resid = " << resid);

         if (resid <= final_norm)
         {
            Update(x, i, H, s, v);
            final_norm = resid;
            final_iter = j;
            converged = true;
            goto finish;
         }

         if (print_options.iterations)
         {
            mfem::out << "   Pass : " << setw(2) << (j-1)/m+1
                     << "   Iteration : " << setw(3) << j
                     << "  ||B r|| = " << resid << '\n';
         }

         Monitor(j, resid, r, x);
      }

      if (print_options.iterations && j <= max_iter)
      {
         mfem::out << "Restarting..." << '\n';
      }

      Update(x, i-1, H, s, v);

      oper->Mult(x, r);
      if (prec)
      {
         subtract(b, r, w);
         prec->Mult(w, r);    // r = M (b - A x)
      }
      else
      {
         subtract(b, r, r);
      }
      beta = Norm(r);         // beta = ||r||
      MFEM_ASSERT(IsFinite(beta), "beta = " << beta);
      if (beta <= final_norm)
      {
         final_norm = beta;
         final_iter = j;
         converged = true;
         goto finish;
      }
   }

   final_norm = beta;
   final_iter = max_iter;
   converged = false;

finish:
   if ((print_options.iterations && converged) || print_options.first_and_last)
   {
      mfem::out << "   Pass : " << setw(2) << (j-1)/m+1
               << "   Iteration : " << setw(3) << final_iter
               << "  ||B r|| = " << final_norm << '\n';
   }
   if (print_options.summary || (print_options.warnings && !converged))
   {
      mfem::out << "GMRES: Number of iterations: " << final_iter << '\n';
   }
   if (print_options.warnings && !converged)
   {
      mfem::out << "GMRES: No convergence!\n";
   }

   Monitor(final_iter, final_norm, r, x, true);

   for (i = 0; i < v.Size(); i++)
   {
      delete v[i];
   }
}