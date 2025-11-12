#pragma once

// Operator to apply action of CV outflow bc
#include "mfem.hpp"
#include "../../Navier-Stokes/ins_bc_containers.hpp"
#include <general/forall.hpp>

using namespace std;

namespace ins_solver_serial
{
class CVBC_Op : public mfem::Operator
{
private:

   // bc arrays
   std::vector<ResistanceBC_T> *resistance_bcs = nullptr;
   bool resistance_bc = false;

   std::vector<WindkesselBC_T> *windkessel_bcs = nullptr;
   bool windkessel_bc = false;

   // data
   mutable bool first_flag = true;

   mfem::FiniteElementSpace *vfes = nullptr;

   // temp
   mutable mfem::Vector temp;

public:
   CVBC_Op(int size, mfem::FiniteElementSpace *vfes_) : mfem::Operator(size), vfes(vfes_) 
   { 
      temp.SetSize(size); // square operator
      temp = 0.0;
   };

   void SetResistanceBC(std::vector<ResistanceBC_T> *resistance_bcs_)
   {
      resistance_bc = true;
      resistance_bcs = resistance_bcs_;
   }

   void SetWindkesselBC(std::vector<WindkesselBC_T> *windkessel_bcs_)
   {
      windkessel_bc = true;
      windkessel_bcs = windkessel_bcs_;
   }

   void Mult(const mfem::Vector &k, mfem::Vector &y) const override;
   void AddMult(const mfem::Vector &x, mfem::Vector &y, const double a=1.0) const override;

   void GetDiag(mfem::Vector &diag, double af);

   void MultGrad(const mfem::Vector &k, mfem::Vector &y) const;
   void AddMultGrad(const mfem::Vector &x, mfem::Vector &y, const double a = 1.0) const;

   virtual ~CVBC_Op() {};
};

}