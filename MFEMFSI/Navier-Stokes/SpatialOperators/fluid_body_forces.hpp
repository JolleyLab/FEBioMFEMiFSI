#pragma once

#include "mfem.hpp"

namespace ins_solver
{

/** Class for domain integration of L(v) := (f, v), where
    f=(f1,...,fn) and v=(v1,...,vn). */
class BodyForceIntegrator : public mfem::LinearFormIntegrator
{
private:
	mfem::Vector shape, Qvec;
	mfem::VectorCoefficient &Q;
	mfem::ConstantCoefficient *coeff;

public:
   /// Constructs a domain integrator with a given VectorCoefficient
   BodyForceIntegrator(mfem::VectorCoefficient &QF)
      : Q(QF), coeff(nullptr) { }
   BodyForceIntegrator(mfem::VectorCoefficient &QF, mfem::ConstantCoefficient &coeff_)
   : Q(QF), coeff(&coeff_) { }

   /** Given a particular Finite Element and a transformation (Tr)
       computes the element right hand side element vector, elvect. */
   virtual void AssembleRHSElementVect(const mfem::FiniteElement &el,
                                       mfem::ElementTransformation &Tr,
                                       mfem::Vector &elvect);

   using mfem::LinearFormIntegrator::AssembleRHSElementVect;
};

}