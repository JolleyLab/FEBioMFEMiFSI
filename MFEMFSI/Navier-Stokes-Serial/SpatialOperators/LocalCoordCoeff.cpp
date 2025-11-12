#include "LocalCoordCoeff.hpp"

using namespace mfem;
using namespace ins_solver_serial;

LocalVectorCoefficient::LocalVectorCoefficient(int vd_) : mfem::VectorCoefficient(vd_)
{
   n_vec.SetSize(vd_);
   t_vec.SetSize(vd_);
   bt_vec.SetSize(vd_);

   n_vec = 0.0;
   t_vec = 0.0;
   bt_vec = 0.0;
};

LocalVectorCoefficient::LocalVectorCoefficient(int vd_, double &vn_, double &v1_, double &v2_) : LocalVectorCoefficient(vd_)
{
   vn = vn_;
   v1 = v1_;
   v2 = v2_;
};

LocalVectorCoefficient::LocalVectorCoefficient(int vd_, double &vn_) : LocalVectorCoefficient(vd_)
{
   vn = vn_;
   only_normal = true;
};

LocalVectorCoefficient::LocalVectorCoefficient(int vd_, double &v1_, double &v2_) : LocalVectorCoefficient(vd_)
{
   v1 = v1_;
   v2 = v2_;
   only_tangential_components = true;
};

LocalVectorCoefficient::LocalVectorCoefficient(int vd_, mfem::Coefficient *vn_coeff_, mfem::Coefficient *v1_coeff_, mfem::Coefficient *v2_coeff_) : LocalVectorCoefficient(vd_)
{
   vn_coeff = vn_coeff_;
   v1_coeff = v1_coeff_;
   v2_coeff = v2_coeff_;
   coeff_flag = true;
};

LocalVectorCoefficient::LocalVectorCoefficient(int vd_, mfem::Coefficient *vn_coeff_) : LocalVectorCoefficient(vd_)
{
   vn_coeff = vn_coeff_;
   coeff_flag = true;
   only_normal = true;
};

LocalVectorCoefficient::LocalVectorCoefficient(int vd_, mfem::Coefficient *v1_coeff_, mfem::Coefficient *v2_coeff_) : LocalVectorCoefficient(vd_)
{
   v1_coeff = v1_coeff_;
   v2_coeff = v2_coeff_;
   coeff_flag = true;
   only_tangential_components = true;
};

void LocalVectorCoefficient::Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip)
{
   V.SetSize(T.GetSpaceDim());
   MFEM_ASSERT(T.GetSpaceDim() > 1, "LocalVectorCoefficient only works for 2D and 3D");
   
   // compute local vectors
   // normal
   CalcOrtho(T.Jacobian(), n_vec);
   double n_vec_mag = n_vec.Norml2();
   n_vec /= n_vec_mag;

   // tangent
   T.Jacobian().GetColumn(0,t_vec);
   double t_vec_mag = t_vec.Norml2();
   t_vec /= t_vec_mag;

   // bi-tangent
   if (T.GetSpaceDim() == 3) 
   {
      T.Jacobian().GetColumn(1,bt_vec);
      double bt_vec_mag = bt_vec.Norml2();
      bt_vec /= bt_vec_mag;
   }

   if (only_normal)
   {
      if(coeff_flag) vn = vn_coeff->Eval(T,ip);

      V.Set(vn,n_vec);
   }
   else if (only_tangential_components)
   {
      MFEM_ABORT("TODO");
      if(coeff_flag) v1 = v1_coeff->Eval(T,ip);
      
      V.Set(v1,t_vec);
      if(T.GetSpaceDim() == 3) 
      {
         if(coeff_flag) v2 = v2_coeff->Eval(T,ip);
         V.Add(v2,bt_vec);
      }
   }
   else
   {
      MFEM_ABORT("TODO");
      if(coeff_flag) vn = vn_coeff->Eval(T,ip);
      if(coeff_flag) v1 = v1_coeff->Eval(T,ip);

      V.Set(vn,n_vec);
      V.Add(v1,t_vec);
      if(T.GetSpaceDim() == 3) 
      {
         if(coeff_flag) v2 = v2_coeff->Eval(T,ip);
         V.Add(v2,bt_vec);
      }
   }
};