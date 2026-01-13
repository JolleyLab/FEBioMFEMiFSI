#pragma once

#include "mfem.hpp"

using VecFuncT = void(const mfem::Vector &x, double t, mfem::Vector &u);
using ScalarFuncT = double(const mfem::Vector &x, double t);
using MatFuncT = void(const mfem::Vector &x, double t, mfem::DenseMatrix &gradu);

/// Container for a displacement boundary condition (dirichlet)
class DisplacementBC_T
{
public:
   DisplacementBC_T(mfem::Array<int> attr, mfem::VectorCoefficient *coeff)
      : attr(attr), coeff(coeff), coeff_array(nullptr)
   {
      all_comp = true;
      component = -1;
   }

   DisplacementBC_T(mfem::Array<int> attr, mfem::Coefficient *coeff_[], int component_)
      : attr(attr), coeff_array(coeff_), coeff(nullptr)
   {
      all_comp = false;
      component = component_;
   }

   DisplacementBC_T(DisplacementBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      this->coeff_array = obj.coeff_array;
      this->all_comp = obj.all_comp;
      this->component = obj.component;
      obj.coeff = nullptr;
      obj.coeff_array = nullptr;
   }

   DisplacementBC_T & operator= (const DisplacementBC_T &obj)
   {
      // move assignment
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      this->coeff_array = obj.coeff_array;
      this->all_comp = obj.all_comp;
      this->component = obj.component;

      return *this;
   };


   ~DisplacementBC_T() { delete coeff; delete coeff_array; }

   mfem::Array<int> attr;
   mfem::VectorCoefficient *coeff;
   mfem::Coefficient **coeff_array;
   bool all_comp;
   int component = -1;
};

/// Container for a traction boundary condition (neumann)
class TractionBC_T
{
public:
   TractionBC_T(mfem::Array<int> attr, mfem::VectorCoefficient *coeff)
      : attr(attr), coeff(coeff), mcoeff(NULL)
   {}

   TractionBC_T(mfem::Array<int> attr, mfem::MatrixCoefficient *mcoeff)
      : attr(attr), mcoeff(mcoeff), coeff(NULL)
   {}

   TractionBC_T(TractionBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      obj.coeff = nullptr;
      this->mcoeff = obj.mcoeff;
      obj.mcoeff = nullptr;
   }

   ~TractionBC_T() 
   { 
      delete coeff;
      delete mcoeff;
   }

   mfem::Array<int> attr;
   mfem::VectorCoefficient *coeff;
   mfem::MatrixCoefficient *mcoeff;
};

/// Container for body force term.
class BodyForce_T
{
public:
   BodyForce_T(mfem::Array<int> attr, mfem::VectorCoefficient *coeff)
      : attr(attr), coeff(coeff)
   {}

   BodyForce_T(BodyForce_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      obj.coeff = nullptr;
   }

   ~BodyForce_T() { delete coeff; }

   mfem::Array<int> attr;
   mfem::VectorCoefficient *coeff;
};

/// Container for a follower pressure loads
class FollowerPressureLoadBC_T
{
public:
   FollowerPressureLoadBC_T(mfem::Array<int> attr, mfem::Coefficient *coeff)
      : attr(attr), coeff(coeff)
   {}

   FollowerPressureLoadBC_T(FollowerPressureLoadBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      obj.coeff = nullptr;
   }
   
   ~FollowerPressureLoadBC_T() { delete coeff; }

   mfem::Array<int> attr;
   mfem::Coefficient *coeff;
};