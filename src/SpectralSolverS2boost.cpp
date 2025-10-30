//
// Created by Peter Zimmerman on 29.10.25.
//

#include <cmath>
#include <cassert>

#include "../include/SpectralSolverS2boost.hpp"
#include "../include/GHPScalars.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fftw3.h>

namespace ublas = boost::numeric::ublas;
using std::vector;

using namespace teuk::literals;


namespace SpecS2boost {

// -----------------------------
// Constructor
// -----------------------------
    SpectralSolver::SpectralSolver(int Nz, int Nphi, bool use_fast_double_fft /*=true*/)
            : Nz_(Nz), Nphi_(Nphi), use_fast_double_fft_(use_fast_double_fft) {
        z_ = legendre_gauss_lobatto(Nz_);
        w_ = barycentric_weights();
        D_ = legendre_diff_matrix(z_); // now returns ublas::matrix<Real>
    }

    // ===========================================================
// dz derivative (ublas multiprecision)
// ===========================================================
    void SpectralSolver::dz(const vector<vector<Complex>> &f,
                            vector<vector<Complex>> &df) const {
        const int Nz = Nz_;
        const int Nphi = Nphi_;

        df.assign(Nz, vector<Complex>(Nphi, Complex(0, 0)));

        // temporary uBLAS vectors
        ublas::vector<Complex> col(Nz), dcol(Nz);

        for (int j = 0; j < Nphi; ++j) {
            // copy column j into uBLAS vector
            for (int i = 0; i < Nz; ++i) col(i) = f[i][j];

            // dcol = D_ * col  (D_ is ublas::matrix<Real> but we stored it as Real -> need Complex multiply)
            // We'll do explicit multiply promoting Real -> Complex
            std::fill(dcol.begin(), dcol.end(), Complex(0, 0));
            for (int i = 0; i < Nz; ++i) {
                for (int k = 0; k < Nz; ++k) {
                    dcol(i) += Complex(D_(i, k)) * col(k);
                }
            }

            // copy back
            for (int i = 0; i < Nz; ++i) df[i][j] = dcol(i);
        }
    } // dz

//   ===================================================
// multiprecision DFT (naive O(N^2)) -- accurate but slow
// used when use_fast_double_fft_ == false
// ===========================================================
    static void dft_mp_line(const vector<Complex> &in, vector<Complex> &out, int sign) {
        // sign = +1 for forward (compute Fourier coeffs), sign = -1 for inverse
        const int N = (int) in.size();
        out.assign(N, Complex(0, 0));
        for (int m = 0; m < N; ++m) {
            Complex sum(0, 0);
            for (int n = 0; n < N; ++n) {
                long double angle = (2.0L * M_PI * m * n) / (long double) N;
                Complex phase = Complex(std::cos(angle), sign * std::sin(angle));
                sum += in[n] * phase;
            }
            if (sign == -1) // inverse DFT normalization
                out[m] = sum / Complex((Real) N);
            else
                out[m] = sum;
        }
    }

// ===========================================================
// dphi derivative: either use FFTW (fast, double precision inside FFT)
// or multiprecision DFT (slow, exact for teuk::Complex)
// ===========================================================
    void SpectralSolver::dphi_fft(const vector<vector<Complex>> &f,
                                  vector<vector<Complex>> &df_dphi) const {
        const int Nz = Nz_;
        const int Nphi = Nphi_;

        df_dphi.assign(Nz, vector<Complex>(Nphi, Complex(0, 0)));

        if (use_fast_double_fft_) {
            // Fast path: cast each phi-line to std::complex<double>, run FFTW,
            // multiply by i*m, inverse FFT, cast result back.
            std::vector<std::complex<double>> tmp_in(Nphi), tmp_out(Nphi);

            for (int iz = 0; iz < Nz; ++iz) {
                // copy & cast to double complex
                for (int j = 0; j < Nphi; ++j) {
                    const Complex &z = f[iz][j];
                    tmp_in[j] = std::complex<double>(static_cast<double>(z.real()),
                                                     static_cast<double>(z.imag()));
                }

                // forward FFT
                fftw_plan p_fwd = fftw_plan_dft_1d(Nphi,
                                                   reinterpret_cast<fftw_complex *>(tmp_in.data()),
                                                   reinterpret_cast<fftw_complex *>(tmp_out.data()),
                                                   FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(p_fwd);
                fftw_destroy_plan(p_fwd);

                // multiply by i*m in Fourier space
                for (int m = 0; m < Nphi; ++m) {
                    int ms = (m <= Nphi / 2) ? m : m - Nphi;
                    tmp_out[m] *= std::complex<double>(0.0, (double) ms);
                }

                // inverse FFT: result back into tmp_in
                fftw_plan p_bwd = fftw_plan_dft_1d(Nphi,
                                                   reinterpret_cast<fftw_complex *>(tmp_out.data()),
                                                   reinterpret_cast<fftw_complex *>(tmp_in.data()),
                                                   FFTW_BACKWARD, FFTW_ESTIMATE);
                fftw_execute(p_bwd);
                fftw_destroy_plan(p_bwd);

                // normalize and cast back to multiprecision Complex
                for (int j = 0; j < Nphi; ++j) {
                    std::complex<double> val = tmp_in[j] / static_cast<double>(Nphi);
                    df_dphi[iz][j] = Complex((Real) val.real(), (Real) val.imag());
                }
            }

        } else {
            // Precise O(N^2) multiprecision DFT path: forward, multiply, inverse (normalized)
            vector<Complex> tmp_freq(Nphi), tmp_ifft(Nphi), tmp_in(Nphi);

            for (int iz = 0; iz < Nz; ++iz) {
                // forward DFT (no 1/N)
                for (int j = 0; j < Nphi; ++j) tmp_in[j] = f[iz][j];

                dft_mp_line(tmp_in, tmp_freq, +1);

                // multiply by i*m (frequency domain)
                for (int m = 0; m < Nphi; ++m) {
                    int ms = (m <= Nphi / 2) ? m : m - Nphi;
                    tmp_freq[m] *= Complex(0, (Real) ms);
                }

                // inverse DFT (normalized inside)
                dft_mp_line(tmp_freq, tmp_ifft, -1);

                // store
                for (int j = 0; j < Nphi; ++j) df_dphi[iz][j] = tmp_ifft[j];
            }
        }
    }

// ===========================================================
// Barycentric weights (multiprecision)
// ===========================================================
    vector<Real> SpectralSolver::barycentric_weights() const {
        int N = (int) z_.size();
        vector<Real> w(N, Real(1));
        for (int j = 0; j < N; ++j) {
            w[j] = Real(1);
            for (int k = 0; k < N; ++k) if (k != j) w[j] /= (z_[j] - z_[k]);
        }
        return w;
    }

// ===========================================================
// Barycentric interpolation & derivative (multiprecision)
// ===========================================================
    std::pair<Complex, Complex> SpectralSolver::barycentric_interp_and_derivative(
            const vector<Complex> &f, const vector<Real> &w, Real z0) const {
        int N = (int) z_.size();
        Complex num(0, 0), num_der(0, 0);
        Real den = Real(0);

        const Real eps = std::numeric_limits<Real>::epsilon() * Real(1e3); // loosened eps for multiprecision
        for (int i = 0; i < N; ++i) {
            Real dz = z0 - z_[i];
            if (std::abs(dz) < eps) { // exact or extremely close to node
                return {f[i],
                        Complex(std::numeric_limits<Real>::quiet_NaN(),
                                std::numeric_limits<Real>::quiet_NaN())};
            }
            Real wi_dz = w[i] / dz;
            num += wi_dz * f[i];
            num_der += wi_dz * (f[i] / dz);
            den += wi_dz;
        }

        Complex f_interp = num / den;
        Complex df_interp = num_der - f_interp * (num_der / num);
        return {f_interp, df_interp};
    }

// ===========================================================
// Legendre polynomial & derivative (multiprecision)
// ===========================================================
    std::pair<Real, Real> SpectralSolver::legendre_P_and_dP(int n, Real x) {
        if (n == 0) return {Real(1), Real(0)};
        if (n == 1) return {x, Real(1)};
        Real Pnm1 = Real(1), Pn = x;
        for (int k = 2; k <= n; ++k) {
            Real Pnp1 = ((Real(2) * k - Real(1)) * x * Pn - (k - 1) * Pnm1) / Real(k);
            Pnm1 = Pn;
            Pn = Pnp1;
        }
        Real dPn = Real(n) / (x * x - Real(1)) * (x * Pn - Pnm1);
        return {Pn, dPn};
    }


// ===========================================================
// Legendre Gauss-Lobatto nodes (multiprecision)
// ===========================================================
    vector<Real> SpectralSolver::legendre_gauss_lobatto(int N){
        vector<Real> x(N);
        x[0] = Real(-1); x[N-1] = Real(1);

        for (int i = 1; i < N-1; ++i) {
            // initial guess using double then casted is fine
            long double xi_ld = -std::cos(M_PI * i / (N - 1));
            Real xi = Real(xi_ld);
            for (int it = 0; it < 40; ++it) {
                auto pr = legendre_P_and_dP(N-1, xi);
                Real P = pr.first;
                Real dP = pr.second;
                Real d2P = (Real(2)*xi*dP - Real((N-1)*N)*P) / (Real(1) - xi*xi);
                Real dx = -dP / d2P;
                xi += dx;
                if (teuk::Fabs(dx) < Real(1e-30)) break;
            }
            x[i] = xi;
        }
        return x;
    }

// ===========================================================
// Legendre differentiation matrix (multiprecision uBLAS matrix)
// ===========================================================
    ublas::matrix<Real> SpectralSolver::legendre_diff_matrix(const vector<Real>& x){
        int N = (int)x.size();
        ublas::matrix<Real> D(N, N);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                D(i,j) = Real(0);

        vector<Real> Pnm1(N);
        for (int i = 0; i < N; ++i) Pnm1[i] = legendre_P_and_dP(N-1, x[i]).first;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i != j) D(i,j) = Pnm1[i] / (Pnm1[j] * (x[i] - x[j]));
            }
        }
        for (int i = 0; i < N; ++i) {
            Real s = Real(0);
            for (int j = 0; j < N; ++j) if (i != j) s += D(i,j);
            D(i,i) = -s;
        }
        return D;
    }


/** ======================================================================
* \function edth
*
* \brief Computes the eth operator (spin-raising derivative) for a GHP scalar.
*
* Applies the spin-weighted derivative along the m-direction of the Newman-Penrose
* tetrad. The returned GHP scalar has its spin weight raised by one:
* (p,q) -> (p+1, q-1).
*
* @param f_in Input GHPScalar on the spectral grid (values stored in 2D vector).
* @return GHPScalar The eth derivative of the input, with updated spin weights.
*/
    GHPField SpectralSolver::edth(const GHPField& f_in) const
    {
        int Nz = Nz_;
        int Nphi = Nphi_;

        // Allocate derivative arrays
        std::vector<std::vector<Complex>> df_dz(Nz, std::vector<Complex>(Nphi)); // polar
        std::vector<std::vector<Complex>> df_dphi(Nz, std::vector<Complex>(Nphi)); // azi

        const auto& f = f_in.values();

        // Compute spectral derivatives
        dz(f, df_dz);         // ∂/∂z
        dphi_fft(f, df_dphi); // ∂/∂φ

        // Compute eth combination
        std::vector<std::vector<Complex>> df_eth(Nz, std::vector<Complex>(Nphi));
        for(int i = 0; i < Nz; ++i) {
            Real z = nodes()[i];
            Real factor = std::sqrt(1.0 - z*z);
            for(int j = 0; j < Nphi; ++j) {
                df_eth[i][j] = -factor * (df_dz[i][j] + Complex(0.0,1.0) * df_dphi[i][j] / factor);
            }
        }

        // Return new field with raised GHP weights
        // Allocate new GHPField of same size
        GHPField edth_f(Nz, Nphi, 0.0, f_in.p() + 1, f_in.q() - 1);

        // Fill in the computed values
        edth_f.set_values(df_eth);
        return edth_f;
    } // edth

} // specS2boost