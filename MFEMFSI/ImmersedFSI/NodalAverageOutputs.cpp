#include "NodalAverageOutputs.h"
#include <mfem.hpp>

namespace nodal_avg
{

void PrincipalValues(const mfem::DenseMatrix &A, double &p1, double &p2, double &p3)
{
   SymmetricEigenvalues(A(0,0), A(1,1), A.Height() == 3 ? A(2,2) : 0.0,
                        A(0,1), A.Height() == 3 ? A(0,2) : 0.0,
                        A.Height() == 3 ? A(1,2) : 0.0, p1, p2, p3);
}

double StressMeasureValue(const mfem::DenseMatrix &stress, StressQuantity measure, int dim)
{
   switch (measure)
   {
      case StressQuantity::VonMises:
      {
         double dev = stress.Trace() * (1.0 / static_cast<double>(dim));
         mfem::DenseMatrix s_dev(stress);
         s_dev(0,0) -= dev;
         s_dev(1,1) -= dev;
         if (dim == 3) s_dev(2,2) -= dev;
         return std::sqrt(1.5 * s_dev.FNorm2());
      }
      case StressQuantity::CauchyXX: return stress(0,0);
      case StressQuantity::CauchyYY: return stress(1,1);
      case StressQuantity::CauchyZZ: return (dim == 3) ? stress(2,2) : 0.0;
      case StressQuantity::CauchyXY: return stress(0,1);
      case StressQuantity::CauchyYZ: return (dim == 3) ? stress(1,2) : 0.0;
      case StressQuantity::CauchyXZ: return (dim == 3) ? stress(0,2) : 0.0;
      case StressQuantity::Principal1:
      case StressQuantity::Principal2:
      case StressQuantity::Principal3:
      {
         double l1,l2,l3;
         PrincipalValues(stress,l1,l2,l3);
         if (measure == StressQuantity::Principal1) return l1;
         if (measure == StressQuantity::Principal2) return l2;
         return l3;
      }
      default: return 0.0;
   }
}

double StrainMeasureValue(const mfem::DenseMatrix &E, StrainQuantity measure, int dim)
{
   switch (measure)
   {
      case StrainQuantity::VonMises:
      {
         double dev = E.Trace() * (1.0 / static_cast<double>(dim));
         mfem::DenseMatrix e_dev(E);
         e_dev(0,0) -= dev;
         e_dev(1,1) -= dev;
         if (dim == 3) e_dev(2,2) -= dev;
         return std::sqrt(1.5 * e_dev.FNorm2());
      }
      case StrainQuantity::Exx: return E(0,0);
      case StrainQuantity::Eyy: return E(1,1);
      case StrainQuantity::Ezz: return (dim == 3) ? E(2,2) : 0.0;
      case StrainQuantity::Exy: return E(0,1);
      case StrainQuantity::Eyz: return (dim == 3) ? E(1,2) : 0.0;
      case StrainQuantity::Exz: return (dim == 3) ? E(0,2) : 0.0;
      case StrainQuantity::Principal1:
      case StrainQuantity::Principal2:
      case StrainQuantity::Principal3:
      {
         double l1,l2,l3;
         PrincipalValues(E,l1,l2,l3);
         if (measure == StrainQuantity::Principal1) return l1;
         if (measure == StrainQuantity::Principal2) return l2;
         return l3;
      }
      default: return 0.0;
   }
}

} // namespace nodal_avg
