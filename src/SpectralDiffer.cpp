//
// Created by Peter Zimmerman on 22.11.25.
//
#include "../include/GHPScalars.hpp"
#include "../include/SpectralGHPField.hpp"
#include "../include/SpectralDiffer.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <cassert>
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
    vector<Real> SpectralDiffer::build_barycentric_weights_z() const {
        int N = (int) z_.size();
        vector<Real> w(N, Real(1));
        for (int j = 0; j < N; ++j) {
            w[j] = Real(1);
            for (int k = 0; k < N; ++k) if (k != j) w[j] /= (z_[j] - z_[k]);
        }
        return w;
    }

    RVector SpectralDiffer::build_barycentric_weights_from_nodes(const RVector &nodes) const {
        auto x = nodes;
        size_t N = x.size();
        RVector w(N, Real(1));
        for (int j = 0; j < N; ++j) {
            w[j] = Real(1);
            for (int k = 0; k < N; ++k) if (k != j) w[j] /= (x[j] - x[k]);
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
            const RVector &x_,          // interpolation nodes
            const CVector &fvec,        // function to interpolate (value and derivative) should be a slice at fixed r or z
            const RVector &w,          // barycentric weights
            Real z0                     // point of interest
    ) {
        int N = (int) x_.size();
        Complex num(0.0, 0.0), num_der(0.0, 0.0);
        Real den = Real(0);

        const Real eps = std::numeric_limits<Real>::epsilon() * Real(1e3);
        for (int i = 0; i < N; ++i) {
            Real dz = z0 - x_[i];
            if (abs(dz) < eps) {
                return {fvec[i],
                        Complex(std::numeric_limits<Real>::quiet_NaN(),
                                std::numeric_limits<Real>::quiet_NaN())};
            }
            Real wi_dz = w[i] / dz;
            num += wi_dz * fvec[i];
            num_der += wi_dz * (fvec[i] / dz);
            den += wi_dz;
        }

        Complex f_interp = num / den;
        Complex df_interp = num_der - f_interp * (num_der / num);
        return {f_interp, df_interp};
    }

/**
* @brief Differentiation using barycentric formula.
*
* df_dz[i] = sum_{j != i} w[j] / w[i] / (z[i]-z[j]) * (f[j]-f[i])
*/
    void SpectralDiffer::dz_barycentric_inplace(const std::span<GHPScalar<Complex>> &f_slice,
                                                std::span<GHPScalar<Complex>> &df_dz,
                                                const std::vector<Real> &w) const {
        if (f_slice.size() != static_cast<size_t>(Nz_) || df_dz.size() != static_cast<size_t>(Nz_))
            throw std::runtime_error("Slice size mismatch in dz_barycentric");

#pragma omp parallel for default(none) shared(f_slice, df_dz, w)
        for (int i = 0; i < Nz_; ++i) {
            Complex sum = Complex(0.0, 0.0);
            for (int j = 0; j < Nz_; ++j) {
                if (i == j) continue;
                sum += w[j] / w[i] * (f_slice[j].value() - f_slice[i].value()) / (z_[i] - z_[j]);
            }
            df_dz[i].value() = sum;
        }
    }

    GHPSpectralField::RSlice SpectralDiffer::dz_barycentric_RSlice(const GHPSpectralField::RSlice &f,
                                                                   GHPSpectralField::RSlice &df_dz,
                                                                   const std::vector<Real> &w) const {
        assert(f.size() == df_dz.size());
        const int N = (int) f.size();
        const Real eps = std::numeric_limits<Real>::epsilon() * Real(1e3);

        for (int i = 0; i < N; ++i) {
            Complex num(0, 0), num_der(0, 0);
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

    void SpectralDiffer::dr_barycentric_inplace(const std::span<GHPScalar<Complex>> &f_slice, // fixed z
                                                std::span<GHPScalar<Complex>> &df_dr,         // output result
                                                const std::vector<Real> &w) const {
        // weights for r

        if (f_slice.size() != static_cast<size_t>(Nr_) || df_dr.size() != static_cast<size_t>(Nr_))
            throw std::runtime_error("Slice size mismatch in dr_barycentric");

#pragma omp parallel for default(none) shared(f_slice, df_dr, w)
        for (int i = 0; i < Nr_; ++i) {
            Complex sum = Complex(0.0, 0.0);
            for (int j = 0; j < Nr_; ++j) {
                if (i == j) continue;
                sum += w[j] / w[i] * (f_slice[j].value() - f_slice[i].value()) / (r_[i] - r_[j]);
            }
            df_dr[i].value() = sum;
        }
    }

    GHPSpectralField::ZSlice SpectralDiffer::dr_barycentric_ZSlice(const GHPSpectralField::ZSlice &f,
                                                                   GHPSpectralField::ZSlice &df_dr,
                                                                   const std::vector<Real> &w) const {
        assert(f.size() == df_dr.size());
        const int N = (int) f.size();
        const Real eps = std::numeric_limits<Real>::epsilon() * Real(1e3);

        for (int i = 0; i < N; ++i) {
            Complex num(0, 0), num_der(0, 0);
            Real den = Real(0);

            Real r0 = r_[i]; // collocation point where derivative is computed

            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                Real dr = r0 - r_[j];
                if (abs(dr) < eps) dr = eps; // avoid divide-by-zero
                Real wi_dr = w[j] / dr;
                num += wi_dr * f[j].value();
                num_der += wi_dr * (f[j].value() / dr);
                den += wi_dr;
            }

            Complex f_val = f[i].value(); // collocation value at i
            Complex df_val = num_der - (num / den) * (num_der / den);

            // wrap into GHPScalar preserving spin weights
            df_dr[i] = GHPScalar<Complex>(df_val, f[i].p(), f[i].q());
        }

        return df_dr;
    }

    vector<Real> SpectralDiffer::build_chebyshev_lobatto_nodes(int N) {
        vector <Real> x(N);

        // Chebyshev–nodes of the second kind:
        // x_j = cos(pi * j / (N - 1))
        for (int j = 0; j < N; ++j)
            x[j] = Real(-cos(M_PI * Real(j) / Real(N - teuk::one)));

        return x;
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
 *   Note that the analytic cancellation in  P'_n(x) = n (x P_n(x) − P_{n−1}(x)) / (x² − 1),
 *   causes loss of significant digits when x→±1, even with multiple precision, bc both the numeratora
 *   and denominartor approach zero as O(1-x). For this reason we implement
 *   an analytic endpoint expansion of the derivative
 *
 *
 * This routine is used heavily in:
 *   - Gauss–Lobatto node construction
 *   - Construction of the Legendre differentiation matrix
 */
    pair<Real, Real> SpectralDiffer::legendre_P_and_dP(unsigned int n, const Real &x) {
        if (n == 0) return {Real(1), Real(0)};
        if (n == 1) return {x, Real(1)};
        Real Pnm1 = Real(1.);    // P_{0}
        Real Pn = x;                // P_{1}

        // Recurrence to compute P_n(x)
        for (int k = 2; k <= n; ++k) {
            Real Pnp1 = ((Real(2) * k - Real(1.)) * x * Pn - Real(k - 1.0) * Pnm1) / Real(k);
            Pnm1 = Pn; // shift up
            Pn = Pnp1; // shift up
        }

        //
        // handle the points x ~ ±1 carefully to avoid loss of precision
        //
        // Multiprecision epsilon
        Real eps = std::numeric_limits<Real>::epsilon();
        // Distance to endpoint scaled by precision
        Real dist = Real(1.0) - abs(x);
        Real tol = sqrt(eps);

        // case 1: away from endpoints
        // safe to evaluate derivative using P'_n(x) = n (x P_n(x) − P_{n−1}(x)) / (x² − 1),
        if (dist > tol) {
            Real denom = x * x - Real(1.);
            Real numer = x * Pn - Pnm1;
            Real dPn = Real(n) * numer / denom;
            return {Pn, dPn};
        }
        // CASE 2: too close to ±1 → use analytic limit
        // P_n'(±1) = ± (-1)*sign*n(n+1)/2
        Real dPn_endpoint = Real(n * (n + 1) / 2.0);

        // (-1)^{n+1} factor
        if (x < 0) {
            if ((n % 2) == 0)
                dPn_endpoint = -dPn_endpoint;
        }

        // P_n(±1)
        Real Pn_endpoint =
                (x >= 0 ? Real(1) : ((n % 2 == 0) ? Real(1) : Real(-1)));

        return {Pn_endpoint, dPn_endpoint};
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
    vector<Real> SpectralDiffer::build_legendre_gauss_lobatto_nodes(int N) {
        vector <Real> x(N); // nodes
        x[0] = Real(-1);
        x[N - 1] = Real(1); // LGL includes endpoints

        auto iMax = 40;
        for (int i = 1; i < N - 1; ++i) {
            Real xi_ld = -cos(M_PI * i / (N - 1)); // splitting of the interval [-1,1] into N pieces
            Real xi = Real(xi_ld);
            for (int it = 0; it < iMax; ++it) {
                auto pr = legendre_P_and_dP(N - 1, xi);
                Real P = pr.first;
                Real dP = pr.second;
                Real d2P = (Real(2) * xi * dP - Real((N - 1) * N) * P) / (Real(1) - xi * xi);
                Real dx = -dP / d2P;
                xi += dx;
                if (abs(dx) < Real(1e-30)) break;
            }
            x[i] = xi;
        }
        return x;
    }

    matrix<Real> SpectralDiffer::build_legendre_diff_matrix(const RVector &x) {
        int N = (int) x.size();
        matrix <Real> D(N, N);
        vector <Real> Pnm1(N);
        for (int i = 0; i < N; ++i) Pnm1[i] = legendre_P_and_dP(N - 1, x[i]).first;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i != j) D(i, j) = Pnm1[i] / (Pnm1[j] * (x[i] - x[j]));
                else D(i, j) = Real(0);
            }
        }
        for (int i = 0; i < N; ++i) {
            Real s = Real(0);
            for (int j = 0; j < N; ++j) if (i != j) s += D(i, j);
            D(i, i) = -s;
        }
        return D;
    }
/**
 * @brief Build Chebyshev differentiation matrix for N nodes x.
 */
    matrix<Real> SpectralDiffer::build_chebyshev_diff_matrix(const RVector &x) {
        int N = static_cast<int>(x.size());
        matrix <Real> D(N, N);
        std::vector<Real> c(N);

        for (int i = 0; i < N; ++i)
            c[i] = (i == 0 || i == N - 1) ? 2.0 : 1.0;
        for (int i = 0; i < N; ++i) c[i] = c[i] * ((i % 2) == 0 ? 1.0 : -1.0);

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i != j)
                    D(i, j) = (c[i] / c[j]) / (x[i] - x[j]);
                else
                    D(i, j) = 0.0;
            }
        }

        for (int i = 0; i < N; ++i) {
            Real s = 0.0;
            for (int j = 0; j < N; ++j) if (i != j) s += D(i, j);
            D(i, i) = -s;
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
 * Legendre–Gauss–Lobatto differentiation matrix Dz_. This performs a dense
 * matrix–vector multiplication:
 *
 *      (df/dz)_i = Σ_j  Dz_(i, j) * f_j ,
 *
 * where each f_j is a GHPScalar and Dz_(i, j) is a real spectral coefficient.
 *
 * The spin-weights (p, q) of the input slice f are preserved in the output
 * by constructing df_dz[i] with the same (p, q). The differentiation matrix
 * itself is spin-neutral, so only the complex scalar components are acted on.
 *
 * Complexity: O(Nz²) for each call, since the full dense differentiation
 * matrix is applied.
 *
 */
    void SpectralDiffer::dz_Dmatrix_RSlice(const GHPSpectralField::RSlice &f,
                                           GHPSpectralField::RSlice &df_dz) const {
        assert(f.size() == df_dz.size());
        for (int i = 0; i < Nz_; ++i) {
            df_dz[i] = GHPScalar<Complex>(teuk::zeroC, f[i].p(), f[i].q());
            for (int j = 0; j < Nz_; ++j) {
                df_dz[i] = df_dz[i] + f[j] * GHPScalar<Complex>(Dz_(i, j), 0, 0);  // matrix multiplication
            }
        }
    }

/**
 * @brief Compute dz using the Legendre differentiation matrix for a single row.
 *
 * @param f_slice Input slice (row of GHPFieldVectorized)
 * @param df_dz  Output slice (preallocated)
 */
    void SpectralDiffer::dz_Dmatrix(const std::span<GHPScalar<Complex>> &f_slice,
                                    std::span<GHPScalar<Complex>> &df_dz) const {

        if (f_slice.size() != static_cast<size_t>(Nz_) || df_dz.size() != static_cast<size_t>(Nz_))
            throw std::runtime_error("Slice size mismatch in dz_Dmatrix");

#pragma omp parallel for default(none) shared(f_slice, df_dz)
        for (int i = 0; i < Nz_; ++i) {
            Complex sum = Complex(0.0, 0.0);
            for (int j = 0; j < Nz_; ++j)
                sum += Dz_(i, j) * f_slice[j].value();
            df_dz[i].value() = sum;
            // Spin weights remain unchanged
        }
    }
/**
 * @brief Compute dr using the Cheb. differentiation matrix for a single row.
 *
 * @param f_slice Input slice (row of GHPFieldVectorized)
 * @param df_dr  Output slice (preallocated)
 */

    void SpectralDiffer::dr_Dmatrix(const std::span<GHPScalar<Complex>> &f_slice,
                                    std::span<GHPScalar<Complex>> &df_dr) const {

        const int N = static_cast<int>(f_slice.size());
        if (df_dr.size() != static_cast<size_t>(N))
            throw std::runtime_error("Slice size mismatch in dr_Dmatrix_slice");

#pragma omp parallel for default(none) shared(f_slice, df_dr, Dr_, N)
        for (int i = 0; i < N; ++i) {
            Complex sum = Complex(0.0, 0.0);
            for (int j = 0; j < N; ++j)
                sum += Dr_(i, j) * f_slice[j].value();
            df_dr[i].value() = sum;
            // spin weights remain unchanged
        }
    }

/**
 * @brief derivative along phi for a single m-mode.
 *
 */
    void SpectralDiffer::dphi_single_m_slice(const std::span<GHPScalar<Complex>> &f_slice,
                                             std::span<GHPScalar<Complex>> &df_dphi_slice,
                                             const Real &m) const {
        if (f_slice.size() != df_dphi_slice.size())
            throw std::runtime_error("Slice size mismatch in dphi_fft_single_m");

        // Example: assume f_slice is already Fourier transformed, multiply by i*m
        int Nz_local = static_cast<int>(f_slice.size());
#pragma omp parallel for default(none) shared(f_slice, df_dphi_slice, m, Nz_local)
        for (int i = 0; i < Nz_local; ++i) {
            df_dphi_slice[i].value() = Complex(0, 1) * f_slice[i].value() * static_cast<Real>(m);
        }
    }

    void SpectralDiffer::dt_single_omega_slice(const std::span<GHPScalar<Complex>> &f_slice,
                                               std::span<GHPScalar<Complex>> &df_dphi_slice,
                                               const Real &omega) const {
        // Same as above; could include w-dependent scaling if needed
        dphi_single_m_slice(f_slice, df_dphi_slice, -omega);
    }


/**

            Real z = lgl_nodes()[i];
            Complex factor = Complex(Sqrt(1.0 - z * z), 0);
            Complex df_dphi = f_RSlice[i].value() * Complex(0, 1) * Real(m);
            out[i].value() = -Complex(1.0 / Sqrt(2.0), 0) * (
                    -factor * df_dz_RSlice[i].value()
                    + Complex(0, 1) * df_dphi / factor
                    - s * z * f_RSlice[i].value() / factor
                    + aw * factor * f_RSlice[i].value());
        */

// -----------------------------------------------------------------------------
// In-place edthH using barycentric dz
// -----------------------------------------------------------------------------
    void SpectralDiffer::edthH_bary_inplace(const GHPSpectralField::RSlice &f,
                                            GHPSpectralField::RSlice &out, Real const& a) const
    {
        // safety
        assert(f.size() == out.size());
        const int N = f.size();
        if (N == 0) return;

        // metadata
        const int p = f[0].p();
        const int q = f[0].q();
        const Complex s = Complex((p - q) * 0.5, 0.0);
        const Real m = Real(f.m());
        const Real aw = a*f.w() * Real(2.*M_PI) / (Nt * dt);

        // Temporary buffer for df/dz (owned here -> safe if out aliases f)
        std::vector<GHPScalar<Complex>> df_dz_buf;
        df_dz_buf.assign(N, GHPScalar<Complex>(Complex(0.0, 0.0), p, q));

        GHPSpectralField::RSlice df_dz_rs(df_dz_buf.data(), N, f.w(), f.m(), f.r());

        // compute derivative in-place into df_dz_buf
        dz_barycentric_RSlice(f, df_dz_rs, wz_); // assumes this signature exists

        // compute edth_H into out (spin weights p+1, q-1)
        const Complex pref = Complex(-1.0 / std::sqrt(2.0), 0.0);

        for (int i = 0; i < N; ++i) {
            const Real z = lgl_nodes()[i];
            const Real factor = std::sqrt(std::max<Real>(0.0, 1.0 - z*z));

            const Complex fval = f[i].value();
            const Complex dfz  = df_dz_rs[i].value();
            const Complex dfphi = Complex(0.0, 1.0) * Real(m) * fval;

            Complex val = pref * (
                    -factor * dfz
                    + Complex(0.0, 1.0) * dfphi / factor
                    - s * z * fval / factor
                    + Complex(aw * factor, 0.0) * fval
            );

            out[i] = GHPScalar<Complex>(val, p + 1, q - 1);
        }
    }

// -----------------------------------------------------------------------------
// In-place edthH using Legendre D-matrix derivative
// -----------------------------------------------------------------------------
    void SpectralDiffer::edthH_dmat_RSlice_inplace(const GHPSpectralField::RSlice &f,
                                            GHPSpectralField::RSlice &out, const Real& a) const
    {
        // safety
        assert(f.size() == out.size());
        const int N = f.size();
        if (N == 0) return;

        // metadata
        const int p = f[0].p();
        const int q = f[0].q();
        const Complex s = Complex((p - q) * 0.5, 0.0); // spin weight
        const Real m = Real(f.m());
        const Real aw = a*f.w() * Real(2.*M_PI) / (Nt * dt);

        // Temporary buffer for df/dz
        std::vector<GHPScalar<Complex>> df_dz_buf;
        df_dz_buf.assign(N, GHPScalar<Complex>(Complex(0.0, 0.0), p, q));

        GHPSpectralField::RSlice df_dz_rs(df_dz_buf.data(), N, f.w(), f.m(), f.r());

        // compute derivative using D-matrix into df_dz_buf
        dz_Dmatrix_RSlice(f, df_dz_rs); // assumes this signature exists

        // compute edth_H into out (spin weights p+1, q-1)
        const Complex pref = Complex(-1.0 / std::sqrt(2.0), 0.0);

        for (int i = 0; i < N; ++i) {
            const Real z = lgl_nodes()[i];
            const Real factor = std::sqrt(std::max<Real>(0.0, 1.0 - z*z));

            const Complex fval = f[i].value();
            const Complex dfz  = df_dz_rs[i].value();
            const Complex dfphi = Complex(0.0, 1.0) * Real(m) * fval;

            Complex val = pref * (
                    -factor * dfz
                    + Complex(0.0, 1.0) * dfphi / factor
                    - s * z * fval / factor
                    + Complex(aw * factor, 0.0) * fval
            );

            out[i] = GHPScalar<Complex>(val, p + 1, q - 1);
        }
    }

    void SpectralDiffer::edthH_inplace_RSliceV(
            const spectral::SpectralGHPVectorized::RSlice &in_RSlice,
            spectral::SpectralGHPVectorized::RSlice &out_RSlice, const Real& aw) const
    {
        // sanity check: sizes must match
        assert(in_RSlice.size() == out_RSlice.size());
        const size_t Nz = in_RSlice.size();

        // extract GHP weights and mode info
        int p = in_RSlice[0].p();
        int q = in_RSlice[0].q();
        Complex s = Complex((p - q) / 2., 0);
        Real m = Real(in_RSlice.m());

        // temporary derivative along z
        std::vector<GHPScalar<Complex>> df_dz(Nz);
        // wrap as RSlice at the same r
        spectral::SpectralGHPVectorized::RSlice df_dz_slice(df_dz.data(), Nz,
                                                            in_RSlice.modes_,
                                                            in_RSlice.r());

        // compute dz derivative
        std::span<GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
        std::span<GHPScalar<Complex>> df_dz_span(df_dz_slice.data_ptr, Nz);

        dz_Dmatrix(in_span, df_dz_span);  // computes derivative along z into df_dz_span

        // main edthH computation
#pragma omp parallel for default(none) shared(Nz,p, q, out_RSlice, in_RSlice, df_dz_slice, m, s, aw)
        for (size_t i = 0; i < Nz; ++i)
        {
            Real z = lgl_nodes()[i];                     // Legendre-Gauss-Lobatto nodes
            Complex factor = Complex(std::sqrt(1.0 - z * z), 0);
            Complex df_dphi = in_RSlice[i].value() * Complex(0, 1) * Real(m);

            out_RSlice[i].value() = -Complex(1.0 / std::sqrt(2.0), 0) * (
                    -factor * df_dz_slice[i].value()
                    + df_dphi / factor * Complex(0, 1)
                    - s * z * in_RSlice[i].value() / factor
                    + aw * factor * in_RSlice[i].value()
            );

            // update GHP weights (edth raises spin by 1)
            out_RSlice[i].set_pq(p + 1, q - 1);
        }
    }

    void SpectralDiffer::edthBarH_inplace_RSliceV(
            const SpectralGHPVectorized::RSlice &in_RSlice,
            SpectralGHPVectorized::RSlice &out_RSlice,
            const Real& aw) const
    {
        // sanity check: sizes must match
        assert(in_RSlice.size() == out_RSlice.size());
        const size_t Nz = in_RSlice.size();

        // extract GHP weights and mode info
        int p = in_RSlice[0].p();
        int q = in_RSlice[0].q();
        Complex s = Complex((p - q) / 2., 0);
        Real m = Real(in_RSlice.m());

        // temporary derivative along z
        std::vector<GHPScalar<Complex>> df_dz(Nz);
        spectral::SpectralGHPVectorized::RSlice df_dz_slice(df_dz.data(), Nz,
                                                            in_RSlice.modes_,
                                                            in_RSlice.r());

        // compute dz derivative
        std::span<GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
        std::span<GHPScalar<Complex>> df_dz_span(df_dz_slice.data_ptr, Nz);
        dz_Dmatrix(in_span, df_dz_span);

        // main edthBar computation
#pragma omp parallel for default(none) shared(Nz,p,q,out_RSlice,in_RSlice,df_dz_slice,m,s,aw)
        for (size_t i = 0; i < Nz; ++i)
        {
            Real z = lgl_nodes()[i];                     // Legendre-Gauss-Lobatto nodes
            Complex factor = Complex(std::sqrt(1.0 - z * z), 0);
            Complex df_dphi = in_RSlice[i].value() * Complex(0, 1) * Real(m);

            // edthBar flips the dz sign and some spin terms compared to edthH
            out_RSlice[i].value() = -Complex(1.0 / std::sqrt(2.0), 0) * (
                   - factor * df_dz_slice[i].value()
                    - df_dphi / factor * Complex(0, 1)
                    + s * z * in_RSlice[i].value() / factor
                    - aw * factor * in_RSlice[i].value()
            );

            // update GHP weights (edthBar lowers spin by 1)
            out_RSlice[i].set_pq(p - 1, q + 1);
        }
    }



} // namespace spectral