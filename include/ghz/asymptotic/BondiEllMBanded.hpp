#ifndef GHZ_NUMERIC_BONDI_ELL_M_BANDED_HPP
#define GHZ_NUMERIC_BONDI_ELL_M_BANDED_HPP

// Small coefficient-space pieces of the Kerr Bondi-gauge reconstruction.
//
// BondiGauge_Adjusted_Generic.nb represents multiplication by cos(theta) on
// spin-weighted spherical harmonics as a real tridiagonal matrix.  The
// adjusted Phi fields then use 1/rho = -r + i a cos(theta).  Keeping this
// operation sparse is important: applying it three times is the only ell
// padding needed by the rho^{-3} term, and it should not turn the reconstruction
// into a dense matrix calculation.

#include "ghz/core/GhzTypes.hpp"

#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace ghz::asymptotic {

using EllMComplex = teuk::Complex;
using EllMReal = teuk::Real;
using EllMVector = Eigen::Matrix<EllMComplex, Eigen::Dynamic, 1>;
using EllMSparseMatrix = Eigen::SparseMatrix<EllMComplex>;

// Coefficients from cos(theta) {}_sY_lm, as used in the notebook:
//   cSM(l,s,m) = sqrt[((l^2-m^2)(l^2-s^2))/(l^2(4l^2-1))]
//   bSM(l,s,m) = -s m/[l(l+1)].
inline EllMReal cos_theta_lower_coefficient(int ell, int spin, int m) {
  if (ell <= 0)
    return EllMReal(0);
  const EllMReal l = EllMReal(ell);
  const EllMReal numerator = (l*l - EllMReal(m*m)) *
                             (l*l - EllMReal(spin*spin));
  const EllMReal denominator = l*l * (EllMReal(4)*l*l - EllMReal(1));
  return numerator > EllMReal(0)
             ? std::sqrt(numerator / denominator)
             : EllMReal(0);
}

inline EllMReal cos_theta_diagonal_coefficient(int ell, int spin, int m) {
  const EllMReal l = EllMReal(ell);
  return -EllMReal(spin*m) / (l * (l + EllMReal(1)));
}

// Rows and columns are ordered according to ells.  The matrix is deliberately
// assembled from triplets rather than through dense temporary storage.
inline EllMSparseMatrix cos_theta_matrix(int spin, int m,
                                         const std::vector<int>& ells) {
  EllMSparseMatrix matrix(static_cast<Eigen::Index>(ells.size()),
                          static_cast<Eigen::Index>(ells.size()));
  std::vector<Eigen::Triplet<EllMComplex>> entries;
  entries.reserve(3 * ells.size());
  std::unordered_map<int, Eigen::Index> index;
  index.reserve(ells.size());
  for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(ells.size()); ++j)
    index.emplace(ells[static_cast<std::size_t>(j)], j);
  for (Eigen::Index col = 0; col < static_cast<Eigen::Index>(ells.size()); ++col) {
    const int lp = ells[static_cast<std::size_t>(col)];
    const auto add = [&](int L, EllMReal value) {
      const auto it = index.find(L);
      if (it != index.end() && value != EllMReal(0))
        entries.emplace_back(it->second, col,
                             EllMComplex(value, EllMReal(0)));
    };
    add(lp + 1, cos_theta_lower_coefficient(lp + 1, spin, m));
    add(lp, cos_theta_diagonal_coefficient(lp, spin, m));
    add(lp - 1, cos_theta_lower_coefficient(lp, spin, m));
  }
  matrix.setFromTriplets(entries.begin(), entries.end());
  return matrix;
}

// Horner evaluation of
//   Phi0 + rho^{-1} Phi1 + rho^{-2} Phi2 + rho^{-3} Phi3,
// with rho^{-1} represented in the spin -2 spherical-harmonic basis by
//   R = -r I + i a C_cos.
//
// The vectors must all have the same ell ordering and length.  The caller
// should include at least three extra ell values when the output is truncated
// at ell_max, as in PhiAdjustmentYVector in the notebook.
inline EllMVector phi_adjustment_horner(const EllMVector& phi0,
                                        const EllMVector& phi1,
                                        const EllMVector& phi2,
                                        const EllMVector& phi3,
                                        EllMReal r, EllMReal a,
                                        const EllMSparseMatrix& cos_matrix) {
  const Eigen::Index n = phi0.size();
  if (phi1.size() != n || phi2.size() != n || phi3.size() != n ||
      cos_matrix.rows() != n || cos_matrix.cols() != n)
    throw std::invalid_argument("Phi vectors and cos(theta) matrix have incompatible sizes");
  EllMSparseMatrix R(n, n);
  R.setIdentity();
  R *= EllMComplex(-r, EllMReal(0));
  R += EllMComplex(0, a) * cos_matrix;
  // This is intentionally written in the same nested order as the notebook;
  // it minimizes sparse products and makes the truncation convention explicit.
  return phi0 + R * phi1 + R * (R * phi2) + R * (R * (R * phi3));
}

} // namespace ghz::asymptotic

#endif
