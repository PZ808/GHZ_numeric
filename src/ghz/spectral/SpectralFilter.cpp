//
// Created by Peter Zimmerman on 04.03.26.
//
#include "ghz/spectral/SpectralFilter.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace spectral {

    // row-major index conventions
    static inline std::size_t idx2(std::size_t i, std::size_t j, std::size_t stride) noexcept {
        return i*stride + j;
    }

    SpectralModalFilter::SpectralModalFilter(std::size_t Nr,
                                             std::size_t Nz,
                                             const std::vector<Real>& z_nodes,
                                             Params params)
            : Nr_(Nr), Nz_(Nz), params_(params)
    {
        if (Nr_ < 2 || Nz_ < 2) throw std::runtime_error("SpectralModalFilter: Nr and Nz must be >= 2");
        if (z_nodes.size() != Nz_) throw std::runtime_error("SpectralModalFilter: z_nodes.size mismatch");

        // ---- Chebyshev cosine table for DCT-I on Lobatto nodes
        cos_r_.resize(Nr_ * Nr_);
        const Real denom_r = Real(Nr_ - 1);
        for (std::size_t k = 0; k < Nr_; ++k) {
            for (std::size_t j = 0; j < Nr_; ++j) {
                cos_r_[idx2(k, j, Nr_)] = std::cos(M_PI * Real(k) * Real(j) / denom_r);
            }
        }
        sigma_x_ = build_sigma(Nr_, params_.alpha_r, params_.p_r);

        // ---- Legendre P_ell(z_i) table
        Pz_.resize(Nz_ * Nz_);
        std::vector<Real> Ptmp(Nz_, 0.0);
        for (std::size_t i = 0; i < Nz_; ++i) {
            legendre_all_(z_nodes[i], Ptmp);
            for (std::size_t ell = 0; ell < Nz_; ++ell) {
                Pz_[idx2(ell, i, Nz_)] = Ptmp[ell];
            }
        }

        // ---- LGL quadrature weights on z-grid
        wz_ = build_lgl_weights_(z_nodes);
        sigma_z_ = build_sigma(Nz_, params_.alpha_z, params_.p_z);
    }

    void SpectralModalFilter::set_params(const Params& params) {
        params_ = params;
        sigma_r_ = build_sigma(Nr_, params_.alpha_r, params_.p_r);
        sigma_z_ = build_sigma(Nz_, params_.alpha_z, params_.p_z);
    }

    std::vector<SpectralModalFilter::Real>
    SpectralModalFilter::build_sigma(std::size_t N, Real alpha, int p)
    {
        std::vector<Real> s(N, Real(1));
        if (N <= 1) return s;
        const Real denom = Real(N - 1);

        // Always preserve low modes; sigma_0 = 1
        s[0] = Real(1);
        for (std::size_t k = 1; k < N; ++k) {
            const Real x = Real(k) / denom;     // in (0,1]
            const Real xp = std::pow(x, Real(p));
            s[k] = std::exp(-alpha * xp);
        }
        return s;
    }

// ------------------------
// Chebyshev (Lobatto) transforms
// ------------------------
//
// Forward: coefficients c_k from nodal values f_j at x_j = cos(pi*j/(N-1)):
//   c_k = (2/(N-1)) [ 0.5 f_0 + sum_{j=1}^{N-2} f_j cos(pi*k*j/(N-1)) + 0.5 f_{N-1} (-1)^k ]
//
// Inverse:
//   f_j = 0.5 c_0 + sum_{k=1}^{N-2} c_k cos(pi*k*j/(N-1)) + 0.5 c_{N-1} (-1)^j
//
    void SpectralModalFilter::cheb_forward_(const std::vector<Complex>& f_j,
                                            std::vector<Complex>& c_k) const
    {
        c_k.assign(Nr_, Complex(0,0));
        const Real scale = Real(2) / Real(Nr_ - 1);

        for (std::size_t k = 0; k < Nr_; ++k) {
            Complex sum = Complex(0,0);
            // endpoints half weight
            sum += Real(0.5) * f_j[0];
            sum += Real(0.5) * f_j[Nr_-1] * ( (k % 2 == 0) ? Real(1) : Real(-1) ); // (-1)^k

            for (std::size_t j = 1; j < Nr_-1; ++j) {
                sum += f_j[j] * cos_x_[idx2(k, j, Nr_)];
            }
            c_k[k] = scale * sum;
        }
    }

    void SpectralModalFilter::cheb_inverse_(const std::vector<Complex>& c_k,
                                            std::vector<Complex>& f_j) const
    {
        f_j.assign(Nr_, Complex(0,0));

        for (std::size_t j = 0; j < Nr_; ++j) {
            Complex sum = Complex(0,0);
            sum += Real(0.5) * c_k[0];
            sum += Real(0.5) * c_k[Nr_-1] * ( (j % 2 == 0) ? Real(1) : Real(-1) ); // (-1)^j

            for (std::size_t k = 1; k < Nr_-1; ++k) {
                sum += c_k[k] * cos_r_[idx2(k, j, Nr_)];
            }
            f_j[j] = sum;
        }
    }

// ------------------------
// Legendre transforms on LGL nodes
// ------------------------
//
// Projection using LGL quadrature weights (exact up to degree 2N-3):
//   a_ell = (2ell+1)/2 * sum_i w_i f_i P_ell(x_i)
//
// Reconstruction:
//   f_i ≈ sum_ell a_ell P_ell(x_i)
//
    void SpectralModalFilter::leg_forward_(const std::vector<Complex>& f_i,
                                           std::vector<Complex>& a_ell) const
    {
        a_ell.assign(Nz_, Complex(0,0));
        for (std::size_t ell = 0; ell < Nz_; ++ell) {
            Complex s = Complex(0,0);
            for (std::size_t i = 0; i < Nz_; ++i) {
                s += wz_[i] * f_i[i] * Pz_[idx2(ell, i, Nz_)];
            }
            const Real norm = Real(0.5) * (Real(2)*Real(ell) + Real(1)); // (2ell+1)/2
            a_ell[ell] = norm * s;
        }
    }

    void SpectralModalFilter::leg_inverse_(const std::vector<Complex>& a_ell,
                                           std::vector<Complex>& f_i) const
    {
        f_i.assign(Nz_, Complex(0,0));
        for (std::size_t i = 0; i < Nz_; ++i) {
            Complex s = Complex(0,0);
            for (std::size_t ell = 0; ell < Nz_; ++ell) {
                s += a_ell[ell] * Pz_[idx2(ell, i, Nz_)];
            }
            f_i[i] = s;
        }
    }

    void SpectralModalFilter::legendre_all_(Real x, std::vector<Real>& P)
    {
        const std::size_t N = P.size();
        if (N == 0) return;
        P[0] = Real(1);
        if (N == 1) return;
        P[1] = x;
        for (std::size_t n = 2; n < N; ++n) {
            const Real nn = Real(n);
            P[n] = ((Real(2)*nn - Real(1))*x*P[n-1] - (nn-Real(1))*P[n-2]) / nn;
        }
    }

    std::vector<SpectralModalFilter::Real>
    SpectralModalFilter::build_lgl_weights_(const std::vector<Real>& x)
    {
        const std::size_t N = x.size();
        if (N < 2) throw std::runtime_error("build_lgl_weights_: N<2");
        std::vector<Real> w(N, Real(0));

        // Need P_{N-1}(x_i)
        std::vector<Real> P(N, 0.0);

        for (std::size_t i = 0; i < N; ++i) {
            legendre_all_(x[i], P);
            const Real PN1 = P[N-1];
            const Real denom = Real(N) * Real(N-1) * (PN1 * PN1);
            w[i] = Real(2) / denom;
        }
        return w;
    }

// ------------------------
// Public filtering methods
// ------------------------

    void SpectralModalFilter::filter_r_inplace(std::vector<Complex>& f) const
    {
        if (f.size() != Nr_ * Nz_) throw std::runtime_error("filter_r_inplace: size mismatch");

        std::vector<Complex> slice(Nr_);
        std::vector<Complex> coeff(Nr_);
        std::vector<Complex> back(Nr_);

        // For each fixed z-index, take r-slice and filter
#pragma omp parallel for if(Nz_>8) default(none) shared(f) firstprivate(slice,coeff,back)
        for (std::ptrdiff_t iz = 0; iz < (std::ptrdiff_t)Nz_; ++iz) {
            // gather
            for (std::size_t ir = 0; ir < Nr_; ++ir) {
                slice[ir] = f[ir * Nz_ + (std::size_t)iz];
            }

            // forward -> filter -> inverse
            cheb_forward_(slice, coeff);
            for (std::size_t k = 0; k < Nr_; ++k) coeff[k] *= sigma_r_[k];
            cheb_inverse_(coeff, back);

            // scatter
            for (std::size_t ir = 0; ir < Nr_; ++ir) {
                f[ir * Nz_ + (std::size_t)iz] = back[ir];
            }
        }
    }

    void SpectralModalFilter::filter_z_inplace(std::vector<Complex>& f) const
    {
        if (f.size() != Nr_ * Nz_) throw std::runtime_error("filter_z_inplace: size mismatch");

        std::vector<Complex> slice(Nz_);
        std::vector<Complex> coeff(Nz_);
        std::vector<Complex> back(Nz_);

        // For each fixed r-index, take z-slice and filter
#pragma omp parallel for if(Nr_>8) default(none) shared(f) firstprivate(slice,coeff,back)
        for (std::ptrdiff_t ir = 0; ir < (std::ptrdiff_t)Nr_; ++ir) {
            // gather
            const std::size_t rbase = (std::size_t)ir * Nz_;
            for (std::size_t iz = 0; iz < Nz_; ++iz) {
                slice[iz] = f[rbase + iz];
            }

            // forward -> filter -> inverse
            leg_forward_(slice, coeff);
            for (std::size_t ell = 0; ell < Nz_; ++ell) coeff[ell] *= sigma_z_[ell];
            leg_inverse_(coeff, back);

            // scatter
            for (std::size_t iz = 0; iz < Nz_; ++iz) {
                f[rbase + iz] = back[iz];
            }
        }
    }

    void SpectralModalFilter::filter_rz_inplace(std::vector<Complex>& f) const
    {
        filter_r_inplace(f);
        filter_z_inplace(f);
    }

} // namespace spectral