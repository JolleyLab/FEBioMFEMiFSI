#pragma once

#include "mfem.hpp"

namespace mfem
{
/// A 0D point finite element
class PointFiniteElementExt : public PointFiniteElement
{
public:
   /// Construct the PointFiniteElement
   PointFiniteElementExt() : PointFiniteElement() {};

   /// Compute the Hessian of second order partial derivatives at @a ip.
   void CalcHessian(const IntegrationPoint &ip,
                    DenseMatrix &h) const override;
};

/// A 1D linear element with nodes on the endpoints
class Linear1DFiniteElementExt : public Linear1DFiniteElement
{
public:
   /// Construct the Linear1DFiniteElement
   Linear1DFiniteElementExt() : Linear1DFiniteElement() {};

   /// Compute the Hessian of second order partial derivatives at @a ip.
   void CalcHessian(const IntegrationPoint &ip,
                    DenseMatrix &h) const override;
};

/// A 2D linear element on triangle with nodes at the vertices of the triangle
class Linear2DFiniteElementExt : public Linear2DFiniteElement
{
public:
   /// Constructor
   Linear2DFiniteElementExt() : Linear2DFiniteElement() {};

   /// Compute the Hessian of second order partial derivatives at @a ip.
   void CalcHessian(const IntegrationPoint &ip,
                    DenseMatrix &h) const override;
};

/// A linear element defined on a triangular prism
class LinearWedgeFiniteElementExt : public LinearWedgeFiniteElement
{
public:
   /// Constructor
   LinearWedgeFiniteElementExt() : LinearWedgeFiniteElement() {};

   /// Compute the Hessian of second order partial derivatives at @a ip.
   void CalcHessian(const IntegrationPoint &ip,
                    DenseMatrix &h) const override;
};

/// A linear element defined on a square pyramid
class LinearPyramidFiniteElementExt : public LinearPyramidFiniteElement
{
public:
   /// Constructor
   LinearPyramidFiniteElementExt() : LinearPyramidFiniteElement() {};

   /// Compute the Hessian of second order partial derivatives at @a ip.
   void CalcHessian(const IntegrationPoint &ip,
                    DenseMatrix &h) const override;
};

/** @brief A 3D linear element on a tetrahedron with nodes at the
    vertices of the tetrahedron */
class Linear3DFiniteElementExt : public Linear3DFiniteElement
{
public:
   /// Constructor
   Linear3DFiniteElementExt() : Linear3DFiniteElement() {};

   /// Compute the Hessian of second order partial derivatives at @a ip.
   void CalcHessian(const IntegrationPoint &ip,
                    DenseMatrix &h) const override;
};

/// A 3D tri-linear element on a cube with nodes at the vertices of the cube
class TriLinear3DFiniteElementExt : public TriLinear3DFiniteElement
{
public:
   /// Constructor
   TriLinear3DFiniteElementExt() : TriLinear3DFiniteElement() {};

   /// Compute the Hessian of second order partial derivatives at @a ip.
   void CalcHessian(const IntegrationPoint &ip,
                    DenseMatrix &h) const override;
};

/// Piecewise-(bi/tri)linear continuous finite elements.
/// with added calc hessian functions
class LinearFECollectionExt : public LinearFECollection
{
private:
   const PointFiniteElementExt PointFE_ext;
   const Linear1DFiniteElementExt SegmentFE_ext;
   const Linear2DFiniteElementExt TriangleFE_ext;
   const BiLinear2DFiniteElement QuadrilateralFE_ext;
   const Linear3DFiniteElementExt TetrahedronFE_ext;
   const TriLinear3DFiniteElementExt ParallelepipedFE_ext;
   const LinearWedgeFiniteElementExt WedgeFE_ext;
   const LinearPyramidFiniteElementExt PyramidFE_ext;
public:
   LinearFECollectionExt() : LinearFECollection() {};

   const FiniteElement *
   FiniteElementForGeometry(Geometry::Type GeomType) const override;
};

} // namespace mfem