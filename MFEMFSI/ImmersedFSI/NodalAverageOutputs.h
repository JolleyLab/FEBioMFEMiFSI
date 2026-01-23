#pragma once

#include <algorithm>
#include <cmath>

// Forward declaration for MFEM types
namespace mfem { class DenseMatrix; }

// MFEM-specific nodal averaging utilities for stress/strain outputs.
// This file is self-contained for the MFEMFSI module.
namespace nodal_avg
{
   enum class StressQuantity
   {
      VonMises,
      CauchyXX,
      CauchyYY,
      CauchyZZ,
      CauchyXY,
      CauchyYZ,
      CauchyXZ,
      Principal1,
      Principal2,
      Principal3
   };

   enum class StrainQuantity
   {
      VonMises,
      Exx,
      Eyy,
      Ezz,
      Exy,
      Eyz,
      Exz,
      Principal1,
      Principal2,
      Principal3
   };

   // Analytic eigenvalues of a symmetric 3x3 matrix (a11,a22,a33,a12,a13,a23).
   // Returns eigenvalues sorted in descending order.
   // see https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
   inline void SymmetricEigenvalues(double a11, double a22, double a33,
                                    double a12, double a13, double a23,
                                    double &l1, double &l2, double &l3)
   {
      const double p1 = a12*a12 + a13*a13 + a23*a23;
      if (p1 == 0.0)
      {
         l1 = a11; l2 = a22; l3 = a33;
      }
      else
      {
         const double trace = a11 + a22 + a33;
         const double q = trace / 3.0;
         const double b11 = a11 - q;
         const double b22 = a22 - q;
         const double b33 = a33 - q;
         const double p2 = b11*b11 + b22*b22 + b33*b33 + 2.0*p1;
         const double p = std::sqrt(p2 / 6.0);

         const double invp = (p > 0.0) ? 1.0/p : 0.0;
         const double b12 = a12 * invp;
         const double b13 = a13 * invp;
         const double b23 = a23 * invp;
         const double c11 = b11 * invp;
         const double c22 = b22 * invp;
         const double c33 = b33 * invp;

         // det(B) for symmetric 3x3
         const double detB = c11*(c22*c33 - b23*b23) - b12*(b12*c33 - b13*b23) + b13*(b12*b23 - c22*b13);
         double r = 0.5 * detB;
         r = std::max(-1.0, std::min(1.0, r));

         const double phi = std::acos(r) / 3.0;
         l1 = q + 2.0*p*std::cos(phi);
         l3 = q + 2.0*p*std::cos(phi + (2.0*M_PI/3.0));
         l2 = trace - l1 - l3;
      }

      // sort descending
      if (l1 < l2) std::swap(l1,l2);
      if (l2 < l3) std::swap(l2,l3);
      if (l1 < l2) std::swap(l1,l2);
   }

   void PrincipalValues(const mfem::DenseMatrix &A, double &p1, double &p2, double &p3);
   double StressMeasureValue(const mfem::DenseMatrix &stress, StressQuantity measure, int dim);
   double StrainMeasureValue(const mfem::DenseMatrix &E, StrainQuantity measure, int dim);
}
