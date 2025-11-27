//
// Created by Peter Zimmerman on 22.11.25.
//

#include <cmath>
#include <cassert>

#include "../include/GHPScalars.hpp"
#include "../include/SpectralGHPField.hpp"
#include "../include/HeldScalars.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <fftw3.h>

#include <limits>

namespace spectral {
    using std::vector;
    using std::pair;


/**
 * @name build_barycentric_weights
 * @return A vector of barycentric weights corresponding to the grid points z_.
 *
 * @brief Compute barycentric interpolation weights for the current node set.
 *
 * Barycentric interpolation expresses the Lagrange interpolant in a numerically
 * stable form that relies on a set of precomputed weights { w_j }. For a given
 * set of interpolation nodes { z_j }, the barycentric weights are defined as:
 *
 *      w_j = 1 / ∏_{k ≠ j} (z_j - z_k)
 *
 * These weights depend only on the nodes, not on the function to be
 * interpolated. Once computed, they can be reused for evaluating the
 * interpolating polynomial and its derivative at arbitrary points.
 *
 */
    vector<Real> SpectralDiffer::build_barycentric_weights() const
    {
        int N = (int)z_.size();
        vector<Real> w(N, Real(1));
        for (int j = 0; j < N; ++j) {
            w[j] = Real(1);
            for (int k = 0; k < N; ++k) if (k != j) w[j] /= (z_[j] - z_[k]);
        }
        return w;
    }

/**
 * @name  barycentric_interp_and_derivative
 * @param fvec  Vector of function values f(z_i) at interpolation nodes.
 * @param w     Barycentric interpolation weights.
 * @param z0    Target point at which to evaluate interpolant and derivative.
 *
 * @return pair ( f(z0), f'(z0) )
 *
 * @brief Performs barycentric Lagrange interpolation and computes the interpolated
 * value and its first derivative at a target point z0.
 *
 * Given:
 *   - fvec[i] : function values f(z_i)
 *   - w[i]    : barycentric weights for the interpolation nodes z_i
 *   - z_[i]   : Gauss–Lobatto (or arbitrary) interpolation nodes
 *
 * The barycentric interpolant is
 *       f(z0) = ( Σ_i w_i f_i / (z0 − z_i) ) / ( Σ_i w_i / (z0 − z_i) ).
 *
 * The derivative is computed using the identity:
 *       f'(z0) =
 *           ( Σ_i w_i f_i / (z0 − z_i)^2 ) / ( Σ_i w_i / (z0 − z_i) )
 *         − f(z0) * ( Σ_i w_i / (z0 − z_i)^2 ) / ( Σ_i w_i / (z0 − z_i) ).
 *
 * Special case:
 *   If z0 coincides with one of the nodes (|z0 − z_i| < eps),
 *   interpolation is exact: f(z0) = f_i and f'(z0) is returned as NaN.
 *
 * This routine is used for:
 *   - Interpolating spectral data to arbitrary points
 *   - Constructing z-derivatives without using the differentiation matrix
 *   - Boundary condition enforcement in spectral codes
.
 */
    pair<Complex, Complex> SpectralDiffer::barycentric_interp_and_derivative(
            const vector<Complex>& fvec,    // function to interpolate (value and derivative)
            const vector<Real>& w,          // barycentric weights
            Real z0                         // point of interest
    ) const
    {
        int N = (int)z_.size();
        Complex num(0.0, 0.0), num_der(0.0, 0.0);
        Real den = Real(0);

        const Real eps = std::numeric_limits<Real>::epsilon() * Real(1e3);
        for (int i = 0; i < N; ++i) {
            Real dz = z0 - z_[i];
            if (abs(dz) < eps) {
                return { fvec[i],
                         Complex(std::numeric_limits<Real>::quiet_NaN(),
                                 std::numeric_limits<Real>::quiet_NaN()) };
            }
            Real wi_dz = w[i] / dz;
            num += wi_dz * fvec[i];
            num_der += wi_dz * (fvec[i] / dz);
            den += wi_dz;
        }

        Complex f_interp = num / den;
        Complex df_interp = num_der - f_interp * (num_der / num);
        return { f_interp, df_interp };
    }

    GHPSpectralField::ZSlice SpectralDiffer::dz_barycentric(
            const GHPSpectralField::ZSlice &f,
            GHPSpectralField::ZSlice &df_dz,
            const std::vector<Real> &w  // barycentric weights
    ) const {
        assert(f.size() == df_dz.size());
        const int N = (int)f.size();
        const Real eps = std::numeric_limits<Real>::epsilon() * Real(1e3);

        for (int i = 0; i < N; ++i) {
            Complex num(0,0), num_der(0,0);
            Real den = Real(0);

            Real z0 = z_[i]; // collocation point where derivative is computed

            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                Real dz = z0 - z_[j];
                if (abs(dz) < eps) dz = eps; // avoid divide-by-zero
                Real wi_dz = w[j] / dz;
                num += wi_dz * f[j].value();
                num_der += wi_dz * (f[j].value() / dz);
                den += wi_dz;
            }

            Complex f_val = f[i].value(); // collocation value at i
            Complex df_val = num_der - (num / den) * (num_der / den);

            // wrap into GHPScalar preserving spin weights
            df_dz[i] = GHPScalar<Complex>(df_val, f[i].p(), f[i].q());
        }

        return df_dz;
    }
/**
 * @name legendre_P_and_dP
 * @param n  Polynomial degree (n ≥ 0).
 * @param x  Evaluation point.
 * @return   A pair (P_n(x), P'_n(x)).
 *
 * @brief Computes the Legendre polynomial P_n(x) and its first derivative P'_n(x)
 * using stable recurrence relations.
 *
 * Definitions:
 *   - P_0(x) = 1
 *   - P_1(x) = x
 *   - Higher-order polynomials satisfy the three-term recurrence:
 *         P_{k+1}(x) =
 *             [(2k + 1) x P_k(x) − k P_{k−1}(x)] / (k + 1)
 *
 * Derivative formula:
 *   The derivative is computed using
 *         P'_n(x) = n (x P_n(x) − P_{n−1}(x)) / (x² − 1),
 *   which is numerically stable except extremely near x = ±1
 *   (these values never occur for interior Gauss–Lobatto nodes).
 *
 * This routine is used heavily in:
 *   - Gauss–Lobatto node construction
 *   - Construction of the Legendre differentiation matrix
 */
    pair<Real, Real> SpectralDiffer::legendre_P_and_dP(int n, Real x)
    {
        if (n == 0) return {Real(1), Real(0)};
        if (n == 1) return {x, Real(1)};
        Real Pnm1 = Real(1), Pn = x;
        for (int k = 2; k <= n; ++k) {
            Real Pnp1 = ((Real(2)*k - Real(1)) * x * Pn - (k - 1) * Pnm1) / Real(k);
            Pnm1 = Pn;
            Pn = Pnp1;
        }
        Real dPn = Real(n) / (x*x - Real(1)) * (x * Pn - Pnm1);
        return { Pn, dPn };
    }

/**
 * @param N  Number of collocation points.
 * @return   A vector containing the N LGL nodes in [-1, 1].
 * @brief Computes the N Legendre–Gauss–Lobatto (LGL) collocation nodes on [-1, 1].
 *
 * LGL nodes include:
 *   - the two boundary points x = -1 and x = +1
 *   - the N−2 interior points satisfying P'_{N-1}(x_i) = 0
 *     (the roots of the derivative of the Legendre polynomial of degree N−1)
 *
 * These nodes are widely used in pseudospectral and collocation schemes
 * because they include the interval boundaries and provide good interpolation
 * properties.
 *
 * Algorithm:
 *   1. The two endpoints are set explicitly.
 *   2. Interior nodes are initialized using Chebyshev points:
 *          x_i ≈ −cos(pi * i / (N−1)).
 *   3. Each interior node is refined via Newton iteration applied to
 *          f(x)  = P'_{N-1}(x)
 *          f'(x) = d²P_{N-1}/dx²
 *      using the identity:
 *          d²P_n/dx² = (2x P'_n − n(n+1) P_n) / (1 − x²).
 *
 * Convergence is typically extremely fast; up to 40 iterations is used as a
 * conservative upper limit.
 */
    vector<Real> SpectralDiffer::build_legendre_gauss_lobatto_nodes(int N)
    {
        vector<Real> x(N); // nodes
        x[0] = Real(-1); x[N-1] = Real(1); // LGL includes endpoints

        auto iMax = 40;
        for (int i = 1; i < N-1; ++i) {
            Real xi_ld = -cos(M_PI*i/(N-1)); // splitting of the interval [-1,1] into N pieces
            Real xi = Real(xi_ld);
            for (int it = 0; it < iMax; ++it) {
                auto pr = legendre_P_and_dP(N-1, xi);
                Real P = pr.first;
                Real dP = pr.second;
                Real d2P = (Real(2)*xi*dP - Real((N-1)*N)*P) / (Real(1) - xi*xi);
                Real dx = -dP / d2P;
                xi += dx;
                if (abs(dx) < Real(1e-30)) break;
            }
            x[i] = xi;
        }
        return x;
    }

// ------------------------------------------------------------
//
// ------------------------------------------------------------
    matrix<Real> SpectralDiffer::build_legendre_diff_matrix(const vector<Real>& x)
    {
        int N = (int)x.size();
        matrix<Real> D(N, N);
        vector<Real> Pnm1(N);
        for (int i = 0; i < N; ++i) Pnm1[i] = legendre_P_and_dP(N-1, x[i]).first;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i != j) D(i,j) = Pnm1[i] / (Pnm1[j] * (x[i] - x[j]));
                else D(i,j) = Real(0);
            }
        }
        for (int i = 0; i < N; ++i) {
            Real s = Real(0);
            for (int j = 0; j < N; ++j) if (i != j) s += D(i,j);
            D(i,i) = -s;
        }
        return D;
    }
/**
 * @name legendre_dz
 * @param f       Input z-slice values of GHPSpectralField values to be differentiated.
 * @param df_dz   Output z-slice values where the derivative will be written.
 *
 * @brief Apply the spectral differentiation matrix in the z-direction to compute d/dz at all grid points of f
 *
 * Computes the z-derivative of a GHP scalar field slice using the precomputed
 * Legendre–Gauss–Lobatto differentiation matrix D_. This performs a dense
 * matrix–vector multiplication:
 *
 *      (df/dz)_i = Σ_j  D_(i, j) * f_j ,
 *
 * where each f_j is a GHPScalar and D_(i, j) is a real spectral coefficient.
 *
 * The spin-weights (p, q) of the input slice f are preserved in the output
 * by constructing df_dz[i] with the same (p, q). The differentiation matrix
 * itself is spin-neutral, so only the complex scalar components are acted on.
 *
 * Complexity: O(Nz²) for each call, since the full dense differentiation
 * matrix is applied.
 *
 */
    void SpectralDiffer::dz_Dmatrix(const GHPSpectralField::ZSlice &f, GHPSpectralField::ZSlice &df_dz) const {
        assert(f.size() == df_dz.size());
        for (int i = 0; i < Nz_; ++i) {
            df_dz[i] = GHPScalar<Complex>(teuk::zeroC, f[i].p(), f[i].q());
            for (int j = 0; j < Nz_; ++j) {
                df_dz[i] = df_dz[i] + f[j] * GHPScalar<Complex>(D_(i,j),0,0);  // matrix multiplication
            }
        }
    }

    // ---------------------------
    // Phi-derivative via FFT for a single ZSlice
    // ---------------------------
    void SpectralDiffer::dphi_fft_single_m(const GHPSpectralField::ZSlice &f_slice,
                                           GHPSpectralField::ZSlice &df_dphi_slice) const
    {
        int m = f_slice.m();  // <-- now it's stored inside
        int Nz = f_slice.size();

        GHPScalar<Complex> im_scalar(Complex(0, m), 0, 0);

        for (int z=0; z<Nz; ++z) df_dphi_slice[z] = f_slice[z] * im_scalar;
    }
    // ---------------------------
    // t-derivative via FFT for a single ZSlice
    // ---------------------------
    void SpectralDiffer::dphi_fft_single_w(const GHPSpectralField::ZSlice &f_slice,
                                           GHPSpectralField::ZSlice &df_dt_slice) const
    {
        int w = f_slice.w();  // <-- now it's stored inside
        int Nz = f_slice.size();

        GHPScalar<Complex> iw_scalar(Complex(0, -w), 0, 0);

        for (int z=0; z<Nz; ++z) df_dt_slice[z] = f_slice[z] * iw_scalar;
    }

    // ---------------------------
    // Edth operator on a single ZSlice
    // ---------------------------
    GHPSpectralField::ZSlice SpectralDiffer::edth(const GHPSpectralField::ZSlice &f) const {

        GHPSpectralField::ZSlice df(f.data_ptr, f.size(), f.w(), f.m(), f.r() );


        // Compute spectral derivatives
        dz_Dmatrix(f, df);               // ∂/∂z using Legendre matrix
        dphi_fft_single_m(f, df); // ∂/∂φ
        dphi_fft_single_w(f, df); // ∂/∂t
        GHPSpectralField::ZSlice df_eth(df.data_ptr, f.size(), f.w(), f.m(), f.r() );


        GHPScalar<Complex> m_scalar(Complex(f.m(), teuk::zero), 0, 0);
        // compute derivative
        for(int i = 0; i < Nz_; ++i) {
            Real z = nodes()[i];
            Real factor = sqrt(1.0 - z*z);
            df_eth[i]= -factor * df[i] - m_scalar * f[i] ;
        }
        return df;
    }



} // namespace spectral
