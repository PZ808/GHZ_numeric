//
// Created by Peter Zimmerman on 04.03.26.
//

#ifndef GHZ_SPECTRAL_SPECTRALFILTER_HPP
#define GHZ_SPECTRAL_SPECTRALFILTER_HPP

#pragma once
#include <vector>
#include <complex>
#include <cstddef>

namespace spectral {

    class SpectralModalFilter {
    public:
        using Real    = double;
        using Complex = std::complex<Real>;

        struct Params {
            // Exponential filter: sigma_k = exp( -alpha * (k/(N-1))^p )
            Real alpha_r = 25.0;
            int  p_r     = 12;

            Real alpha_z = 25.0;
            int  p_z     = 12;

            // If true, also damps k=0..?? No: we keep sigma_0=1 always.
            // (kept for future extension)
        };

        // z_nodes must be LGL nodes in [-1,1] (your Legendre-Gauss-Lobatto grid)
        explicit SpectralModalFilter(std::size_t Nr,
                                     std::size_t Nz,
                                     const std::vector<Real>& z_nodes,
                                     Params params = {});

        void set_params(const Params& params);
        [[nodiscard]] const Params& params() const { return params_; }

        // Filter a 2D field stored row-major: idx = ir*Nz + iz.
        // These apply filtering in-place along the chosen dimension.
        void filter_r_inplace(std::vector<Complex>& f) const; // for each iz, filter over ir
        void filter_z_inplace(std::vector<Complex>& f) const; // for each ir, filter over iz

        // Convenience: filter in both dimensions (r then z)
        void filter_rz_inplace(std::vector<Complex>& f) const;

    private:
        std::size_t Nr_{0}, Nz_{0};
        Params params_;

        // --- Chebyshev (x) precomputations (DCT-I tables) ---
        // cos_x_[k*Nr + j] = cos(pi * k * j / (Nr-1)), k,j=0..Nr-1
        std::vector<Real> cos_x_;
        std::vector<Real> sigma_x_;

        // --- Legendre (z) precomputations ---
        // Pz_[ell*Nz + i] = P_ell(z_i), ell,i=0..Nz-1
        std::vector<Real> Pz_;
        std::vector<Real> wz_;      // LGL quadrature weights
        std::vector<Real> sigma_z_;

        static std::vector<Real> build_sigma(std::size_t N, Real alpha, int p);

        // Chebyshev transforms for Lobatto nodes (DCT-I-like)
        void cheb_forward_(const std::vector<Complex>& f_j, std::vector<Complex>& c_k) const;
        void cheb_inverse_(const std::vector<Complex>& c_k, std::vector<Complex>& f_j) const;

        // Legendre LGL projection and reconstruction
        void leg_forward_(const std::vector<Complex>& f_i, std::vector<Complex>& a_ell) const;
        void leg_inverse_(const std::vector<Complex>& a_ell, std::vector<Complex>& f_i) const;

        // Legendre helper: compute P_n(z) for all n up to N-1 (stable recurrence)
        static void legendre_all_(Real x, std::vector<Real>& P);

        // LGL weights: w_i = 2/(N(N-1) P_{N-1}(z_i)^2)
        static std::vector<Real> build_lgl_weights_(const std::vector<Real>& x);
    };

} // namespace spectral

#endif //GHZ_SPECTRAL_SPECTRALFILTER_HPP
