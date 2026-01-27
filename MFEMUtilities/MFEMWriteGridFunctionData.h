#pragma once

#include <FECore/FEElement.h>
#include <FECore/FEMaterialPoint.h>
#include <functional>
#include <vector>
#include <cmath>

#include "mfem.hpp"
#include "../MFEMImSolid/MFEMNodalAverageOutputs.h"

template <class T> void computeNodalProjectedValues(FEElement& e, mfem::GridFunction& gf, mfem::Array<int>& zones_per_vdof, mfem::Array<int>& vdofs, std::function<T(const FEMaterialPoint&)> var)
{
	T si[FEElement::MAX_INTPOINTS];
	T sn[FEElement::MAX_NODES];

   int ne = e.Nodes();
   int ni = e.GaussPoints();
   int ndofs = vdofs.Size();

   // get the integration point values
   for (int k = 0; k<ni; ++k)
   {
      FEMaterialPoint& mp = *e.GetMaterialPoint(k);
      T s = var(mp);
      si[k] = s;
   }

   // project to nodes
   e.project_to_nodes(si, sn);

   // add data to gf and count zones
   for (int j=0; j<ndofs; ++j)
   {
      gf(vdofs[j]) += sn[j]; 
      zones_per_vdof[vdofs[j]]++;
   }     
}

template <class T>
void computeNodalProjectedTensorMeasureValues(
   FEElement& e,
   mfem::GridFunction& gf,
   mfem::Array<int>& zones_per_vdof,
   mfem::Array<int>& vdofs,
   std::function<mat3ds(const FEMaterialPoint&)> tensor_at_ip,
   std::function<double(const mat3ds&)> eval_measure)
{
   const int ne = e.Nodes();
   const int ni = e.GaussPoints();
   const int ndofs = vdofs.Size();

   // Gather tensor components at integration points
   double si_xx[FEElement::MAX_INTPOINTS];
   double si_yy[FEElement::MAX_INTPOINTS];
   double si_zz[FEElement::MAX_INTPOINTS];
   double si_xy[FEElement::MAX_INTPOINTS];
   double si_yz[FEElement::MAX_INTPOINTS];
   double si_xz[FEElement::MAX_INTPOINTS];

   for (int k = 0; k < ni; ++k)
   {
      FEMaterialPoint& mp = *e.GetMaterialPoint(k);
      const mat3ds S = tensor_at_ip(mp);
      si_xx[k] = S.xx();
      si_yy[k] = S.yy();
      si_zz[k] = S.zz();
      si_xy[k] = S.xy();
      si_yz[k] = S.yz();
      si_xz[k] = S.xz();
   }

   // Project each component to nodes
   double sn_xx[FEElement::MAX_NODES];
   double sn_yy[FEElement::MAX_NODES];
   double sn_zz[FEElement::MAX_NODES];
   double sn_xy[FEElement::MAX_NODES];
   double sn_yz[FEElement::MAX_NODES];
   double sn_xz[FEElement::MAX_NODES];

   e.project_to_nodes(si_xx, sn_xx);
   e.project_to_nodes(si_yy, sn_yy);
   e.project_to_nodes(si_zz, sn_zz);
   e.project_to_nodes(si_xy, sn_xy);
   e.project_to_nodes(si_yz, sn_yz);
   e.project_to_nodes(si_xz, sn_xz);

   // Compute desired measure at nodes and accumulate
   for (int j = 0; j < ndofs; ++j)
   {
      const int node = j;
      mat3ds snode(sn_xx[node], sn_yy[node], sn_zz[node], sn_xy[node], sn_yz[node], sn_xz[node]);
      double val = eval_measure(snode);
      gf(vdofs[j]) += (T)val;
      zones_per_vdof[vdofs[j]]++;
   }
}

// generate all exponent tuples for monomials of total degree <= order in dim variables
inline void generate_exponents(int dim, int order, int idx, int remaining, std::vector<int>& current, std::vector<std::vector<int>>& exps)
{
   if (idx == dim - 1)
   {
      current[idx] = remaining;
      exps.push_back(current);
      return;
   }
   for (int k = 0; k <= remaining; ++k)
   {
      current[idx] = k;
      generate_exponents(dim, order, idx + 1, remaining - k, current, exps);
   }
}

// evaluate a monomial at xi given exponent tuple
inline double eval_monomial(const mfem::Vector& xi, const std::vector<int>& exp)
{
   double val = 1.0;
   for (int d = 0; d < xi.Size(); ++d)
   {
      int p = exp[d];
      if (p == 0) continue;
      double base = xi[d];
      // fast integer power
      double powv = 1.0;
      for (int i = 0; i < p; ++i) powv *= base;
      val *= powv;
   }
   return val;
}

// solve linear system M x = b via Gaussian elimination with partial pivoting
inline bool solve_linear_system(std::vector<std::vector<double>>& M, std::vector<double>& b, std::vector<double>& x)
{
   const int n = (int)M.size();
   if (n == 0 || (int)M[0].size() != n || (int)b.size() != n) return false;
   // Augmented matrix
   for (int i = 0; i < n; ++i) M[i].push_back(b[i]);
   // Elimination
   for (int col = 0; col < n; ++col)
   {
      // Pivot
      int piv = col;
      double best = std::fabs(M[piv][col]);
      for (int r = col + 1; r < n; ++r)
      {
         double v = std::fabs(M[r][col]);
         if (v > best)
         {
            best = v; piv = r;
         }
      }
      if (best == 0.0) return false; // singular
      if (piv != col) std::swap(M[piv], M[col]);
      // Normalize row
      double diag = M[col][col];
      for (int c = col; c <= n; ++c) M[col][c] /= diag;
      // Eliminate below
      for (int r = col + 1; r < n; ++r)
      {
         double factor = M[r][col];
         if (factor == 0.0) continue;
         for (int c = col; c <= n; ++c) M[r][c] -= factor * M[col][c];
      }
   }
   // Back substitution
   x.assign(n, 0.0);
   for (int i = n - 1; i >= 0; --i)
   {
      double sum = M[i][n];
      for (int c = i + 1; c < n; ++c) sum -= M[i][c] * x[c];
      x[i] = sum; // M[i][i] is 1 after normalization
   }
   return true;
}

// Fit a polynomial (in local element coordinates) of given order to IP values, then evaluate at nodes
template <class T>
void computeNodalPolynomialFittedValues(
   FEElement& e,
   const mfem::FiniteElement& mfem_e,
   mfem::GridFunction& gf,
   mfem::Array<int>& zones_per_vdof,
   mfem::Array<int>& vdofs,
   int poly_order,
   std::function<T(const FEMaterialPoint&)> var)
{
   const int ne = e.Nodes();
   const int ni = e.GaussPoints();
   const int ndofs = vdofs.Size();
   const int dim = gf.FESpace()->GetMesh()->SpaceDimension();
   const mfem::IntegrationRule *ir = &mfem_e.GetNodes();

   // Collect IP values y(k) and coordinates
   std::vector<mfem::Vector> xi_ips;
   xi_ips.reserve(ni);
   std::vector<double> y(ni);

   FEElementTraits* traits = e.GetTraits();
   FESolidElementTraits* solid_traits = dynamic_cast<FESolidElementTraits*>(traits);
   if (!solid_traits) MFEM_ABORT("computeNodalPolynomialFittedValues: only solid elements are supported.");
   for (int k = 0; k < ni; ++k)
   {
      mfem::Vector xi(dim);
      xi(0) = solid_traits->gr[k];
      xi(1) = solid_traits->gs[k];
      if (dim == 3) xi(2) = solid_traits->gt[k];
      xi_ips.push_back(xi);
      FEMaterialPoint& mp = *e.GetMaterialPoint(k);
      y[k] = (double)var(mp);
   }

   // Build monomial exponent set up to poly_order
   std::vector<std::vector<int>> exps;
   for (int deg = 0; deg <= poly_order; ++deg)
   {
      std::vector<int> cur(dim, 0);
      generate_exponents(dim, deg, 0, deg, cur, exps);
   }
   const int nterms = (int)exps.size();

   // Assemble normal equations: ATA * c = ATy
   std::vector<std::vector<double>> ATA(nterms, std::vector<double>(nterms, 0.0));
   std::vector<double> ATy(nterms, 0.0);
   for (int k = 0; k < ni; ++k)
   {
      // phi vector at xi_k
      std::vector<double> phi(nterms);
      for (int t = 0; t < nterms; ++t) phi[t] = eval_monomial(xi_ips[k], exps[t]);
      // Update ATA and ATy
      for (int i = 0; i < nterms; ++i)
      {
         ATy[i] += phi[i] * y[k];
         const double pi = phi[i];
         for (int j = i; j < nterms; ++j)
         {
            ATA[i][j] += pi * phi[j];
         }
      }
   }
   // Symmetrize lower triangle
   for (int i = 0; i < nterms; ++i)
      for (int j = 0; j < i; ++j)
         ATA[i][j] = ATA[j][i];

   // Solve for coefficients c
   std::vector<double> c;
   bool ok = solve_linear_system(ATA, ATy, c);

   T sn[FEElement::MAX_NODES];
   if (!ok)
   {
      MFEM_ABORT("computeNodalPolynomialFittedValues: failed to solve for polynomial coefficients.");
      // // Fallback: use existing node projection
      // T si[FEElement::MAX_INTPOINTS];
      // for (int k = 0; k < ni; ++k) si[k] = (T)y[k];
      // e.project_to_nodes(si, sn);
   }
   else
   {
      // Evaluate polynomial at node local coordinates
      for (int i = 0; i < ne; ++i)
      {
         mfem::Vector xi(dim);
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         xi(0) = ip.x;
         xi(1) = ip.y;
         if (dim == 3) xi(2) = ip.z;
         double val = 0.0;
         for (int t = 0; t < nterms; ++t) val += c[t] * eval_monomial(xi, exps[t]);
         sn[i] = (T)val;
      }
   }

   // Accumulate into GridFunction and zone counts
   for (int j = 0; j < ndofs; ++j)
   {
      gf(vdofs[j]) += sn[j];
      zones_per_vdof[vdofs[j]]++;
   }
}

template <class T>
void computeNodalPolynomialFittedTensorMeasureValues(
   FEElement& e,
   const mfem::FiniteElement& mfem_e,
   mfem::GridFunction& gf,
   mfem::Array<int>& zones_per_vdof,
   mfem::Array<int>& vdofs,
   int poly_order,
   std::function<mat3ds(const FEMaterialPoint&)> tensor_at_ip,
   std::function<double(const mat3ds&)> eval_measure)
{
   const int ne = e.Nodes();
   const int ni = e.GaussPoints();
   const int ndofs = vdofs.Size();
   const int dim = gf.FESpace()->GetMesh()->SpaceDimension();
   const mfem::IntegrationRule* ir = &mfem_e.GetNodes();

   // Collect IP local coords (reference space)
   std::vector<mfem::Vector> xi_ips;
   xi_ips.reserve(ni);

   FEElementTraits* traits = e.GetTraits();
   FESolidElementTraits* solid_traits = dynamic_cast<FESolidElementTraits*>(traits);
   if (!solid_traits) MFEM_ABORT("computeNodalPolynomialFittedTensorMeasureValues: only solid elements are supported.");
   for (int k = 0; k < ni; ++k)
   {
      mfem::Vector xi(dim);
      xi(0) = solid_traits->gr[k];
      xi(1) = solid_traits->gs[k];
      if (dim == 3) xi(2) = solid_traits->gt[k];
      xi_ips.push_back(xi);
   }

   // Build monomial exponent set up to poly_order
   std::vector<std::vector<int>> exps;
   for (int deg = 0; deg <= poly_order; ++deg)
   {
      std::vector<int> cur(dim, 0);
      generate_exponents(dim, deg, 0, deg, cur, exps);
   }
   const int nterms = (int)exps.size();

   // Precompute phi(xi_k) and assemble ATA once
   std::vector<std::vector<double>> phi_all(ni, std::vector<double>(nterms));
   std::vector<std::vector<double>> ATA(nterms, std::vector<double>(nterms, 0.0));
   for (int k = 0; k < ni; ++k)
   {
      for (int t = 0; t < nterms; ++t)
      {
         phi_all[k][t] = eval_monomial(xi_ips[k], exps[t]);
      }
      // Update ATA using symmetry
      for (int i = 0; i < nterms; ++i)
      {
         const double pi = phi_all[k][i];
         for (int j = i; j < nterms; ++j)
         {
            ATA[i][j] += pi * phi_all[k][j];
         }
      }
   }
   for (int i = 0; i < nterms; ++i)
      for (int j = 0; j < i; ++j)
         ATA[i][j] = ATA[j][i];

   // Gather IP tensor components
   std::vector<double> y_xx(ni), y_yy(ni), y_zz(ni), y_xy(ni), y_yz(ni), y_xz(ni);
   for (int k = 0; k < ni; ++k)
   {
      FEMaterialPoint& mp = *e.GetMaterialPoint(k);
      const mat3ds S = tensor_at_ip(mp);
      y_xx[k] = S.xx();
      y_yy[k] = S.yy();
      y_zz[k] = S.zz();
      y_xy[k] = S.xy();
      y_yz[k] = S.yz();
      y_xz[k] = S.xz();
   }

   // Build ATy per component
   auto build_ATy = [&](const std::vector<double>& y) {
      std::vector<double> ATy(nterms, 0.0);
      for (int k = 0; k < ni; ++k)
      {
         const double yk = y[k];
         for (int i = 0; i < nterms; ++i) ATy[i] += phi_all[k][i] * yk;
      }
      return ATy;
   };

   std::vector<double> ATy_xx = build_ATy(y_xx);
   std::vector<double> ATy_yy = build_ATy(y_yy);
   std::vector<double> ATy_zz = build_ATy(y_zz);
   std::vector<double> ATy_xy = build_ATy(y_xy);
   std::vector<double> ATy_yz = build_ATy(y_yz);
   std::vector<double> ATy_xz = build_ATy(y_xz);

   // Solve for coefficients per component
   std::vector<double> c_xx, c_yy, c_zz, c_xy, c_yz, c_xz;
   {
      auto M = ATA; if (!solve_linear_system(M, ATy_xx, c_xx)) MFEM_ABORT("Polynomial fit failed (xx)");
   }
   {
      auto M = ATA; if (!solve_linear_system(M, ATy_yy, c_yy)) MFEM_ABORT("Polynomial fit failed (yy)");
   }
   {
      auto M = ATA; if (!solve_linear_system(M, ATy_zz, c_zz)) MFEM_ABORT("Polynomial fit failed (zz)");
   }
   {
      auto M = ATA; if (!solve_linear_system(M, ATy_xy, c_xy)) MFEM_ABORT("Polynomial fit failed (xy)");
   }
   {
      auto M = ATA; if (!solve_linear_system(M, ATy_yz, c_yz)) MFEM_ABORT("Polynomial fit failed (yz)");
   }
   {
      auto M = ATA; if (!solve_linear_system(M, ATy_xz, c_xz)) MFEM_ABORT("Polynomial fit failed (xz)");
   }

   // Evaluate fitted polynomials at node local coordinates
   double sn_xx[FEElement::MAX_NODES];
   double sn_yy[FEElement::MAX_NODES];
   double sn_zz[FEElement::MAX_NODES];
   double sn_xy[FEElement::MAX_NODES];
   double sn_yz[FEElement::MAX_NODES];
   double sn_xz[FEElement::MAX_NODES];

   for (int i = 0; i < ne; ++i)
   {
      mfem::Vector xi(dim);
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      xi(0) = ip.x; xi(1) = ip.y; if (dim == 3) xi(2) = ip.z;
      // Evaluate phi at node
      double phi_node_sum_xx = 0.0, phi_node_sum_yy = 0.0, phi_node_sum_zz = 0.0;
      double phi_node_sum_xy = 0.0, phi_node_sum_yz = 0.0, phi_node_sum_xz = 0.0;
      for (int t = 0; t < nterms; ++t)
      {
         const double phi_t = eval_monomial(xi, exps[t]);
         phi_node_sum_xx += c_xx[t] * phi_t;
         phi_node_sum_yy += c_yy[t] * phi_t;
         phi_node_sum_zz += c_zz[t] * phi_t;
         phi_node_sum_xy += c_xy[t] * phi_t;
         phi_node_sum_yz += c_yz[t] * phi_t;
         phi_node_sum_xz += c_xz[t] * phi_t;
      }
      sn_xx[i] = phi_node_sum_xx;
      sn_yy[i] = phi_node_sum_yy;
      sn_zz[i] = phi_node_sum_zz;
      sn_xy[i] = phi_node_sum_xy;
      sn_yz[i] = phi_node_sum_yz;
      sn_xz[i] = phi_node_sum_xz;
   }

   // Compute desired measure at nodes and accumulate
   for (int j = 0; j < ndofs; ++j)
   {
      const int node = j;
      mat3ds snode(sn_xx[node], sn_yy[node], sn_zz[node], sn_xy[node], sn_yz[node], sn_xz[node]);
      double val = eval_measure(snode);
      gf(vdofs[j]) += (T)val;
      zones_per_vdof[vdofs[j]]++;
   }
}