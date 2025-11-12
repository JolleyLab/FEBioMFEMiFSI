#pragma once

#include "mfem.hpp"

namespace ins_solver_serial
{

// Vector coefficient that sets the value of a vector based on user provided the normal, tangent, and binormal components
class LocalVectorCoefficient : public mfem::VectorCoefficient
{   
private:
   bool only_normal = false;
   bool only_tangential_components = false;

   double vn = 0.0; // normal
   double v1 = 0.0; // tangent
   double v2 = 0.0; // bi-tangent

   // coefficients for each component
   mfem::Coefficient *vn_coeff = nullptr;
   mfem::Coefficient *v1_coeff = nullptr;
   mfem::Coefficient *v2_coeff = nullptr;
   bool coeff_flag = false;

   // vectors
   mfem::Vector n_vec;
   mfem::Vector t_vec;
   mfem::Vector bt_vec;

public:
   LocalVectorCoefficient(int vd_);
   LocalVectorCoefficient(int vd_, double &vn_, double &v1_, double &v2_);
   LocalVectorCoefficient(int vd_, double &vn_);
   LocalVectorCoefficient(int vd_, double &v1_, double &v2_);
   LocalVectorCoefficient(int vd_, mfem::Coefficient *vn_coeff_, mfem::Coefficient *v1_coeff_, mfem::Coefficient *v2_coeff_);
   LocalVectorCoefficient(int vd_, mfem::Coefficient *vn_coeff_);
   LocalVectorCoefficient(int vd_, mfem::Coefficient *v1_coeff_, mfem::Coefficient *v2_coeff_);

   void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
                     const mfem::IntegrationPoint &ip) override;

   virtual ~LocalVectorCoefficient() {};
};

} // ins_solver_serial