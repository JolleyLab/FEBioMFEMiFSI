// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "linear_element_ext.hpp"
#include <config/config.hpp>

using namespace mfem;

void PointFiniteElementExt::CalcHessian(const IntegrationPoint &ip, DenseMatrix &h) const
{
   // nothing to compute
}

void Linear1DFiniteElementExt::CalcHessian(const IntegrationPoint &ip, DenseMatrix &h) const
{
   h(0,0) = 0.0;
   h(1,0) = 0.0;
}

void Linear2DFiniteElementExt::CalcHessian(const IntegrationPoint &ip, DenseMatrix &h) const
{
   h = 0.0;
}

void LinearWedgeFiniteElementExt::CalcHessian(const IntegrationPoint &ip, DenseMatrix &h) const
{
   // CHECK THIS
   h(0,0) = 0.;   h(0,1) =  0.;   h(0,2) =  1.;  h(0,3) =  0.;  h(0,4) =  1.;  h(0,5) =  0.;
   h(1,0) = 0.;   h(1,1) =  0.;   h(1,2) = -1.;  h(1,3) =  0.;  h(1,4) =  0.;  h(1,5) =  0.;
   h(2,0) = 0.;   h(2,1) =  0.;   h(2,2) =  0.;  h(2,3) =  0.;  h(2,4) = -1.;  h(2,5) =  0.;
   h(3,0) = 0.;   h(3,1) =  0.;   h(3,2) = -1.;  h(3,3) =  0.;  h(3,4) = -1.;  h(3,5) =  0.;
   h(4,0) = 0.;   h(4,1) =  0.;   h(4,2) =  1.;  h(4,3) =  0.;  h(4,4) =  0.;  h(4,5) =  0.;
   h(5,0) = 0.;   h(5,1) =  0.;   h(5,2) =  0.;  h(5,3) =  0.;  h(5,4) =  1.;  h(5,5) =  0.;
}

void LinearPyramidFiniteElementExt::CalcHessian(const IntegrationPoint &ip, DenseMatrix &h) const
{
   MFEM_ABORT("Not implemented");
}

void Linear3DFiniteElementExt::CalcHessian(const IntegrationPoint &ip, DenseMatrix &h) const
{
   h = 0.0;
}

void TriLinear3DFiniteElementExt::CalcHessian(const IntegrationPoint &ip, DenseMatrix &h) const
{
   const double x = ip.x, y = ip.y, z = ip.z;
   const double ox = 1.-x, oy = 1.-y, oz = 1.-z;

   h(0,0) = 0.;   h(0,1) =  oz;   h(0,2) =  oy;  h(0,3) =  0.;  h(0,4) =  ox;  h(0,5) =  0.;
   h(1,0) = 0.;   h(1,1) = -oz;   h(1,2) = -oy;  h(1,3) =  0.;  h(1,4) =   x;  h(1,5) =  0.;
   h(2,0) = 0.;   h(2,1) =  oz;   h(2,2) =  -y;  h(2,3) =  0.;  h(2,4) =  -x;  h(2,5) =  0.;
   h(3,0) = 0.;   h(3,1) = -oz;   h(3,2) =   y;  h(3,3) =  0.;  h(3,4) = -ox;  h(3,5) =  0.;
   h(4,0) = 0.;   h(4,1) =   z;   h(4,2) = -oy;  h(4,3) =  0.;  h(4,4) = -ox;  h(4,5) =  0.;
   h(5,0) = 0.;   h(5,1) =  -z;   h(5,2) =  oy;  h(5,3) =  0.;  h(5,4) =  -x;  h(5,5) =  0.;
   h(6,0) = 0.;   h(6,1) =   z;   h(6,2) =   y;  h(6,3) =  0.;  h(6,4) =   x;  h(6,5) =  0.;
   h(7,0) = 0.;   h(7,1) =  -z;   h(7,2) =  -y;  h(7,3) =  0.;  h(7,4) =  ox;  h(7,5) =  0.;
}

const FiniteElement *
LinearFECollectionExt::FiniteElementForGeometry(Geometry::Type GeomType) const
{
   switch (GeomType)
   {
      case Geometry::POINT:       return &PointFE_ext;
      case Geometry::SEGMENT:     return &SegmentFE_ext;
      case Geometry::TRIANGLE:    return &TriangleFE_ext;
      case Geometry::SQUARE:      return &QuadrilateralFE_ext;
      case Geometry::TETRAHEDRON: return &TetrahedronFE_ext;
      case Geometry::CUBE:        return &ParallelepipedFE_ext;
      case Geometry::PRISM:       return &WedgeFE_ext;
      case Geometry::PYRAMID:     return &PyramidFE_ext;
      default:
         if (error_mode == RETURN_NULL) { return nullptr; }
         mfem_error ("LinearFECollectionExt: unknown geometry type.");
   }
   return &SegmentFE_ext; // Make some compilers happy
}