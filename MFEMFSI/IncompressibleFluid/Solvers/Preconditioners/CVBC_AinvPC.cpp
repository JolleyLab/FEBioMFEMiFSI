#include "CVBC_AinvPC.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

void CVBC_AinvPC::SetOperator(const Operator &op)
{
   // to satisfy base class requirement
}

// Operator application
void CVBC_AinvPC::Mult (const Vector & x,Vector & y) const
{   
   // Ainv
   // y.Set(1.0,x);
   // y /= *A3d_diag;
   {
      const bool use_dev = false;
      auto y_vec = y.ReadWrite(use_dev);
      auto x_vec = x.Read(use_dev);
      auto ad_vec = A3d_diag->Read(use_dev);
      const int N = x.Size();
      mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[i] = x_vec[i] / ad_vec[i]; });
   }

   // CVBC contributions (the linearized contributions!)
   if (resistance_bc)
   {
      for (auto &rbc : *resistance_bcs)
      {
         
         Vector *s = rbc.Get_S_vec();

         // Step 1: compute inner product
         // double q_loc = Ad_s * x;
         double q_loc = 0.0;
         for (int i=0; i<rbc.bdr_true_dofs.Size(); i++)
         {
            q_loc += (*s)(i) * x(rbc.bdr_true_dofs[i]) / (*A3d_diag)(rbc.bdr_true_dofs[i]); 
         }
         #ifndef MFEM_USE_MPI
         double q_glob = q_loc;
         #else
         double q_glob = 0.0;
         MPI_Allreduce(&q_loc,&q_glob,1,MPITypeMap<real_t>::mpi_type,MPI_SUM,dynamic_cast<ParFiniteElementSpace*>(vfes)->GetComm());
         #endif

         // Step 2: compute 1 + diag_factor*||Ad_inv_s||^2

         // compute ||Ad_inv_s||^2
         // code from Vector norml2();
         auto data = s->GetData();
         auto data_ad = A3d_diag->GetData();
         double scale = 0.0;
         double sum = 0.0;
         for (int i = 0; i < s->Size(); i++)
         {
            if ((data[i] / data_ad[rbc.bdr_true_dofs[i]]) != 0.0)
            {
               const double absdata = fabs(data[i] / data_ad[rbc.bdr_true_dofs[i]]);
               if (scale <= absdata)
               {
                  const double sqr_arg = scale / absdata;
                  sum = 1.0 + sum * (sqr_arg * sqr_arg);
                  scale = absdata;
                  continue;
               } // end if scale <= absdata
               const double sqr_arg = absdata / scale;
               sum += (sqr_arg * sqr_arg); // else scale > absdata
            } // end if data[i] != 0
         }

         // RB 1-16-24 (need global norm, not local norm!)
         #ifndef MFE_USE_MPI
         double norm_Ad_s =  scale * sqrt(sum);
         double norm_Ad_s_glob = norm_Ad_s;
         #else
         double norm_Ad_s =  scale * sqrt(sum);
         double norm_Ad_s_glob = 0.0;
         MPI_Allreduce(&norm_Ad_s,&norm_Ad_s_glob,1,MPITypeMap<real_t>::mpi_type,MPI_SUM,vfes->GetComm());
         #endif

         // compute the factor
         double fac = 1.0 + rbc.diag_factor*norm_Ad_s_glob;

         // Step 3: scale output vector
         // y.Add(-rbc.diag_factor * q_glob / fac, Ad_s);
         // code from Vector Add()
         const double fac1 = -rbc.diag_factor * q_glob / fac;
         const bool use_dev = false;
         auto y_vec = y.ReadWrite(use_dev);
         auto s_vec = s->Read(use_dev);
         auto ad_vec = A3d_diag->Read(use_dev);
         const int N = rbc.bdr_true_dofs.Size();
         const auto idx = rbc.bdr_true_dofs.Read();
         mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += fac1 * s_vec[i] / ad_vec[idx[i]]; });
      }
   }

   if (windkessel_bc)
   {
      for (auto &wbc : *windkessel_bcs)
      {
         
         Vector *s = wbc.Get_S_vec();

         // Step 1: compute inner product
         // double q_loc = Ad_s * x;
         double q_loc = 0.0;
         for (int i=0; i<wbc.bdr_true_dofs.Size(); i++)
         {
            q_loc += (*s)(i) * x(wbc.bdr_true_dofs[i]) / (*A3d_diag)(wbc.bdr_true_dofs[i]); 
         }
         #ifndef MFEM_USE_MPI
         double q_glob = q_loc;
         #else
         double q_glob = 0.0;
         MPI_Allreduce(&q_loc,&q_glob,1,MPITypeMap<real_t>::mpi_type,MPI_SUM,dynamic_cast<ParFiniteElementSpace*>(vfes)->GetComm());
         #endif

         // Step 2: compute 1 + diag_factor*||Ad_inv_s||^2

         // compute ||Ad_inv_s||^2
         // code from Vector norml2();
         auto data = s->GetData();
         auto data_ad = A3d_diag->GetData();
         double scale = 0.0;
         double sum = 0.0;
         for (int i = 0; i < s->Size(); i++)
         {
            if ((data[i] / data_ad[wbc.bdr_true_dofs[i]]) != 0.0)
            {
               const double absdata = fabs(data[i] / data_ad[wbc.bdr_true_dofs[i]]);
               if (scale <= absdata)
               {
                  const double sqr_arg = scale / absdata;
                  sum = 1.0 + sum * (sqr_arg * sqr_arg);
                  scale = absdata;
                  continue;
               } // end if scale <= absdata
               const double sqr_arg = absdata / scale;
               sum += (sqr_arg * sqr_arg); // else scale > absdata
            } // end if data[i] != 0
         }
         // RB 1-16-24 (need global norm, not local norm!)
         #ifndef MFE_USE_MPI
         double norm_Ad_s =  scale * sqrt(sum);
         double norm_Ad_s_glob = norm_Ad_s;
         #else
         double norm_Ad_s =  scale * sqrt(sum);
         double norm_Ad_s_glob = 0.0;
         MPI_Allreduce(&norm_Ad_s,&norm_Ad_s_glob,1,MPITypeMap<real_t>::mpi_type,MPI_SUM,vfes->GetComm());
         #endif

         // compute the factor
         double fac = 1.0 + wbc.diag_factor*norm_Ad_s_glob;

         // Step 3: scale output vector
         // y.Add(-wbc.diag_factor * q_glob / fac, Ad_s);
         // code from Vector Add()
         const double fac1 = -wbc.diag_factor * q_glob / fac;
         const int N = wbc.bdr_true_dofs.Size();
         const auto idx = wbc.bdr_true_dofs.Read();
         const bool use_dev = false;
         auto y_vec = y.ReadWrite(use_dev);
         auto s_vec = s->Read(use_dev);
         auto ad_vec = A3d_diag->Read(use_dev);
         mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE (int i) { y_vec[idx[i]] += fac1 * s_vec[i] / ad_vec[idx[i]]; });
      }
   }      
}

void CVBC_AinvPC::MultTranspose(const Vector & x,Vector & y) const
{
   MFEM_ABORT("Error: To Do");
}