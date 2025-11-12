#include "CV_Op.hpp"

using namespace std;
using namespace mfem;
using namespace ins_solver_serial;


void CVBC_Op::Mult(const Vector &k, Vector &y) const
{
   MFEM_ABORT("not implemented");
}

// applies (v,n) * f(Q_a) * (u,n)
void CVBC_Op::AddMult(const Vector &x, Vector &y, const double a) const
{
   // resistance bcs apply
   if (resistance_bc)
   {
      for (auto &rbc : *resistance_bcs)
      {
         
         Vector *s = rbc.Get_S_vec();

         // Step 1: compute flow rate
         // only using bdr dofs
         double q_loc = 0.0;
         for (int i=0; i<rbc.bdr_true_dofs.Size(); i++)
         {
            q_loc += x(rbc.bdr_true_dofs[i]) * (*s)(i);
         }
         #ifndef MFEM_USE_MPI
         double q_glob = q_loc;
         #else
         double q_glob = 0.0;
         MPI_Allreduce(&q_loc,&q_glob,1,MPITypeMap<real_t>::mpi_type,MPI_SUM,dynamic_cast<ParFiniteElementSpace*>(vfes)->GetComm());
         #endif

         // Step 2: scale output vector
         const int N = rbc.bdr_true_dofs.Size();
         const auto idx = rbc.bdr_true_dofs.Read();
         auto y_vec = y.ReadWrite();
         auto s_vec = s->Read();
         const double fac = (rbc.R * q_glob + rbc.p0) * a;
         bool use_dev = false;
         mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += fac * s_vec[i]; });
      }
   }

   if (windkessel_bc)
   {
      for (auto &wbc : *windkessel_bcs)
      {      
         Vector *s = wbc.Get_S_vec();

         // exp(-t_n+af / tau) * (- R*Qa(0) - Pd(0)) + Pd(t_n+af)
         // Note: exp(-t_n+af / tau) * P(0) is applied later on if P(0) is nonzero
         double fac0 = 0.0;
         if (wbc.p_d)
         {  
            fac0 = (-wbc.R*wbc.init_flow_rate - wbc.p_d(0))*exp(-(wbc.t_npalphaf)/(wbc.Rd*wbc.C)) + wbc.p_d(wbc.t_npalphaf);
         }
         else
         {
            fac0 = (-wbc.R*wbc.init_flow_rate - wbc.p_d_val)*exp(-(wbc.t_npalphaf)/(wbc.Rd*wbc.C)) + wbc.p_d_val;
         }         

         // hist terms
         double hist_terms = wbc.integ_hist_terms;

         // current term
         // compute current flow rate
         // only using bdr dofs
         double q_loc = 0.0;
         for (int i=0; i<wbc.bdr_true_dofs.Size(); i++)
         {
            q_loc += x(wbc.bdr_true_dofs[i]) * (*s)(i);
         }
         #ifndef MFEM_USE_MPI
         double curr_flow_rate = q_loc;
         #else
         double curr_flow_rate = 0.0;
         MPI_Allreduce(&q_loc,&curr_flow_rate,1,MPITypeMap<real_t>::mpi_type,MPI_SUM,dynamic_cast<ParFiniteElementSpace*>(vfes)->GetComm());
         #endif

         // fac0 = exp(-t_n+af / tau) * (- R*Qa(0) - Pd(0)) + Pd(t_n+af) + R*Q_npaf
         fac0 += wbc.R*curr_flow_rate;

         // midpoint rule approx
         double t_mid = 0.5*(wbc.t_npalphaf + wbc.t_n);
         double t_diff = wbc.t_npalphaf - wbc.t_n;

         wbc.integ_curr_terms = t_diff * (exp(-(wbc.t_npalphaf - t_mid)/(wbc.Rd*wbc.C)) / wbc.C ) * (((t_mid - wbc.t_n)/(wbc.t_npalphaf - wbc.t_n)) * curr_flow_rate 
                              + ((t_mid - wbc.t_npalphaf)/(wbc.t_n - wbc.t_npalphaf)) * wbc.q_n);
         
         // update flow rate at t_n+alpha_f
         wbc.q_npalphaf = curr_flow_rate;

         double curr_term = wbc.integ_curr_terms;

         // set diagonal factor for f'(Q_a) for use later on
         wbc.diag_factor = wbc.R + t_diff * (exp(-(wbc.t_npalphaf - t_mid)/(wbc.Rd * wbc.C)) / wbc.C ) * ((t_mid - wbc.t_n)/(wbc.t_npalphaf - wbc.t_n));

         // Scale output vector
         // y.Add(fac0 + hist_terms + curr_term, *s);
         const double fac1 = fac0 + hist_terms + curr_term;
         const double a_fac = a;
         const int N = wbc.bdr_true_dofs.Size();
         const auto idx = wbc.bdr_true_dofs.Read();
         const bool use_dev = false;
         auto y_vec = y.ReadWrite(use_dev);
         auto s_vec = s->Read(use_dev);
         if (wbc.p_flag)
         {
            // exp(-t_n+af / tau) * P(0) * s_vec
            const double time_fac = exp(-(wbc.t_npalphaf)/(wbc.Rd*wbc.C));
            Vector *p_vec = wbc.Get_p_vec();
            auto p_vec_temp = p_vec->Read(use_dev);
            mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += a_fac * fac1 * s_vec[i] + a_fac * time_fac * p_vec_temp[i]; });
         }
         else
         {
            mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += a_fac * fac1 * s_vec[i]; });
         }
      }
   }
}

// gets diagonal of tangent stiffness contribution (Abc)
// f'(Q_a) * si * si * af
// af = time integration factor
void CVBC_Op::GetDiag(Vector &diag, double af)
{
   first_flag = true;

   diag = 0.0; // zero initially
   if (resistance_bc)
   {
      // TO DO: if not time dependent only form once!
      for (auto &rbc : *resistance_bcs)
      {
         
         Vector *s = rbc.Get_S_vec();

         rbc.diag_factor = rbc.R * af;

         const bool use_dev = false;
         auto y_vec = diag.ReadWrite(use_dev);
         auto s_vec = s->Read(use_dev);
         const double Rval = rbc.diag_factor;
         const int N = rbc.bdr_true_dofs.Size();
         const auto idx = rbc.bdr_true_dofs.Read();
         if (first_flag)
         {
            mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] = Rval * s_vec[i] * s_vec[i]; });
            first_flag = false;
         }
         else
         {
            mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += Rval * s_vec[i] * s_vec[i]; });
         }
      }
   }

   if (windkessel_bc)
   {
      for (auto &wbc : *windkessel_bcs)
      {
         
         Vector *s = wbc.Get_S_vec();

         wbc.diag_factor *= af;

         const int N = wbc.bdr_true_dofs.Size();
         const auto idx = wbc.bdr_true_dofs.Read();
         const bool use_dev = false;
         auto y_vec = diag.ReadWrite(use_dev);
         auto s_vec = s->Read(use_dev);
         const double fval = wbc.diag_factor;
         if (first_flag)
         {
            mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] = fval * s_vec[i] * s_vec[i]; });
            first_flag = false;
         }
         else
         {
            mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += fval * s_vec[i] * s_vec[i]; });
         }
      }
   }

   first_flag = true;
}

void CVBC_Op::MultGrad(const Vector &k, Vector &y) const
{
   MFEM_ABORT("not implemented");
}

// applies (v,n) * f'(Q_a) * (u,n)
void CVBC_Op::AddMultGrad(const Vector &x, Vector &y, const double a) const
{
   // resistance bcs apply
   if (resistance_bc)
   {
      for (auto &rbc : *resistance_bcs)
      {
         Vector *s = rbc.Get_S_vec();

         // Step 1: compute flow rate
         // only using bdr dofs
         double q_loc = 0.0;
         for (int i=0; i<rbc.bdr_true_dofs.Size(); i++)
         {
            q_loc += x(rbc.bdr_true_dofs[i]) * (*s)(i);
         }
         #ifndef MFEM_USE_MPI
         double q_glob = q_loc;
         #else
         double q_glob = 0.0;
         MPI_Allreduce(&q_loc,&q_glob,1,MPITypeMap<real_t>::mpi_type,MPI_SUM,dynamic_cast<ParFiniteElementSpace*>(vfes)->GetComm());
         #endif

         // scale output vector
         const int N = rbc.bdr_true_dofs.Size();
         const auto idx = rbc.bdr_true_dofs.Read();
         auto y_vec = y.ReadWrite();
         auto s_vec = s->Read();
         const double fac = rbc.R * q_glob * a;
         bool use_dev = false;
         mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += fac * s_vec[i]; });
      }
   }

   if (windkessel_bc)
   {
      for (auto &wbc : *windkessel_bcs)
      {      
         Vector *s = wbc.Get_S_vec();

         // compute current flow rate
         // only using bdr dofs
         double q_loc = 0.0;
         for (int i=0; i<wbc.bdr_true_dofs.Size(); i++)
         {
            q_loc += x(wbc.bdr_true_dofs[i]) * (*s)(i);
         }
         #ifndef MFEM_USE_MPI
         double curr_flow_rate = q_loc;
         #else
         double curr_flow_rate = 0.0;
         MPI_Allreduce(&q_loc,&curr_flow_rate,1,MPITypeMap<real_t>::mpi_type,MPI_SUM,dynamic_cast<ParFiniteElementSpace*>(vfes)->GetComm());
         #endif

         // f'(Q_a)
         double fac0 = wbc.diag_factor;

         // Scale output vector
         const int N = wbc.bdr_true_dofs.Size();
         const auto idx = wbc.bdr_true_dofs.Read();
         auto y_vec = y.ReadWrite();
         auto s_vec = s->Read();
         const double fac = fac0*curr_flow_rate * a;
         bool use_dev = false;
         mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += fac * s_vec[i]; });
      }
   }
}