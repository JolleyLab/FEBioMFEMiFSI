#pragma once

#include "mfem.hpp"

namespace ins_solver
{

// surface tractions
// Class for boundary integration of the linear form g = < hv, v >
class TractionIntegrator : public mfem::LinearFormIntegrator
{
   mfem::Vector shape;
   mfem::VectorCoefficient *Q = nullptr;
   mfem::MatrixCoefficient *MQ = nullptr;
   mfem::Coefficient *coeff = nullptr;
   mfem::Coefficient *scalar_coeff = nullptr;
   bool pressure = false; // if true, scalar_coeff represents a pressure value
   int oa, ob;
public:
   /// Constructs a boundary integrator with a given Coefficient QG
   TractionIntegrator(mfem::VectorCoefficient &QG, int a = 1, int b = 1)
      : Q(&QG), MQ(NULL), oa(a), ob(b), coeff(nullptr) { };
   TractionIntegrator(mfem::MatrixCoefficient &MQG, int a = 1, int b = 1)
      : Q(NULL), MQ(&MQG), oa(a), ob(b), coeff(nullptr) { };
   TractionIntegrator(mfem::VectorCoefficient &QG, mfem::Coefficient &coeff_, int a = 1, int b = 1)
      : Q(&QG), MQ(NULL), oa(a), ob(b), coeff(&coeff_) { };
   TractionIntegrator(mfem::MatrixCoefficient &MQG, mfem::Coefficient &coeff_, int a = 1, int b = 1)
      : Q(NULL), MQ(&MQG), oa(a), ob(b), coeff(&coeff_) { };
   TractionIntegrator(mfem::Coefficient &scalar_, int a = 1, int b = 1)
      : Q(NULL), MQ(NULL), oa(a), ob(b), coeff(NULL), scalar_coeff(&scalar_) { };
   TractionIntegrator(mfem::Coefficient &scalar_, bool pressure_, int a = 1, int b = 1)
      : Q(NULL), MQ(NULL), oa(a), ob(b), coeff(NULL), scalar_coeff(&scalar_), pressure(pressure_) { };

   virtual void AssembleRHSElementVect(const mfem::FiniteElement &el,
                                       mfem::ElementTransformation &Tr,
                                       mfem::Vector &elvect);

   using mfem::LinearFormIntegrator::AssembleRHSElementVect;
};

}