//
// Created by Peter Zimmerman on 28.11.25.
//

#ifndef GHZ_NUMERIC_SPECTRALDIFFER_HPP
#define GHZ_NUMERIC_SPECTRALDIFFER_HPP

#include "ghz/core/GhzTypes.hpp"
#include "ghz/core/VectorsGHZ.hpp"
#include "SpectralGHPFieldVectorized.hpp"

#include <vector>
#include <complex>
#include <functional>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <complex>
#include <type_traits>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace spectral {

    using Complex = teuk::Complex;
    using Real = teuk::Real;
    using boost::numeric::ublas::matrix;
    using RVector = std::vector<Real>;
    using CVector = std::vector<Complex>;
    using CMat = std::vector<std::vector<Complex>>;

    /**
     * @class SpectralDiffer
     * @brief Spectral differentiation utilities for functions defined on
     * Legendre–Gauss–Lobatto (LGL) and Chebyshev–Lobatto (CL) nodes.
     *
     * This class provides methods to compute spectral differentiation matrices
     * and perform differentiation operations on functions sampled at LGL and CL nodes.
     *
     */
    class SpectralDiffer {

    public:
        // Constructor: build nodes, weights, differentiation matrix
        SpectralDiffer(size_t Nz, size_t Nr) : Nz_(Nz), Nr_(Nr) {
            z_ = build_legendre_gauss_lobatto_nodes(Nz_); // generate nodes z[i]
            r_ = build_chebyshev_lobatto_nodes(Nr_); // generate nodes r[i]
            wz_ = build_barycentric_weights_LGL(z_);
            wr_ = build_barycentric_weights_ChebLobatto_from_nodes(r_);
            Dz_ = build_legendre_diff_matrix(z_);
            Dr_ = build_chebyshev_diff_matrix(r_);
        }

        // Getters
        [[nodiscard]] const size_t Nr() const { return Nr_; }
        [[nodiscard]] const size_t Nz() const { return Nz_; }
        [[nodiscard]] const RVector &lgl_nodes() const { return z_; }
        [[nodiscard]] const RVector &cl_nodes() const { return r_; }
        [[nodiscard]] const RVector &z_weights() const { return wz_; }
        [[nodiscard]] const RVector &r_weights() const { return wr_; }

        [[nodiscard]] const matrix<Real> &Dz_matrix() const { return Dz_; }
        [[nodiscard]] const matrix<Real> &Dr_matrix() const { return Dr_; }

        // Barycentric utilities (declare; implement elsewhere or below)

        [[nodiscard]] RVector build_generic_barycentric_weights_from_nodes(const RVector &nodes) const;
        [[nodiscard]] RVector build_barycentric_weights_ChebLobatto_from_nodes(const RVector &x) const;
        static RVector build_barycentric_weights_LGL(const RVector& lgl_nodes);

        [[nodiscard]] static std::pair<Complex, Complex> barycentric_interp_and_derivative(const RVector &nodes,
                                                                                           const CVector &f,
                                                                                           const RVector &w,
                                                                                           Real x0);

        static std::pair<Real, Real> legendre_P_and_dP(unsigned int n, const Real &x);

        static std::vector<Real> build_legendre_gauss_lobatto_nodes(int N);

        static std::vector<Real> build_chebyshev_lobatto_nodes(int N);

        static matrix<Real> build_legendre_diff_matrix(const RVector &x);

        static matrix<Real> build_chebyshev_diff_matrix(const RVector &x);


        // dz on a single r=const slice (in-place into df_dz)
        void dz_Dmatrix_RSlice(const SpectralGHPVectorized::RSlice &f_RSlice,
                               SpectralGHPVectorized::RSlice &dfdz_RSlice_into) const;

        void dr_Dmatrix_ZSlice(const spectral::SpectralGHPVectorized::ZSlice &f_ZSlice,
                               SpectralGHPVectorized::ZSlice &df_dr) const;

        void dz_Dmatrix(std::span<const GHPScalar<Complex>> f_span_rconst,
                        std::span<GHPScalar<Complex>> df_dz) const;

        void dr_Dmatrix(std::span<const GHPScalar<Complex>> f_span_zconst,
                        std::span<GHPScalar<Complex>> df_dr) const;

        // main kernels
        void dz_barycentric_inplace(std::span<const GHPScalar<Complex>> f_span_const,
                                    std::span<GHPScalar<Complex>> df) const;

        void dr_barycentric_inplace(std::span<const GHPScalar<Complex>> f_span_const,
                                    std::span<GHPScalar<Complex>> df) const;

        // wrappers
        void dz_barycentric_RSlice(const SpectralGHPVectorized::RSlice &f_RSlice,
                                   SpectralGHPVectorized::RSlice &dfdz_RSlice) const;
        void dr_barycentric_ZSlice(const SpectralGHPVectorized::ZSlice &f,
                                                            SpectralGHPVectorized::ZSlice &df_dr) const;


        static void dphi_single_m_slice(std::span<const GHPScalar<Complex>> f_slice,
                                        std::span<GHPScalar<Complex>> df_dphi_slice,
                                        const Real &m);

        void dt_single_omega_slice(std::span<const GHPScalar<Complex>> f_slice,
                                   std::span<GHPScalar<Complex>> df_dt_slice,
                                   const Real &omega) const;

        // fft routines
        //void dphi_fft(const CMat& f, CMat& df_dphi) const;
        //void dphi_fft_direct(const CMat&  f_hat, CMat&  df_dphi_hat) const;

    private:
        size_t Nz_, Nr_;                  // Number of nodes
        std::vector<Real> z_;            // LGL nodes [-1,1]
        std::vector<Real> r_;            // CL nodes [-1,1]
        std::vector<Real> wz_, wr_;  // barycentric weights
        matrix<Real> Dz_;                // Legendre differentiation matrix
        matrix<Real> Dr_;                // Chebyshev differentiation matrix
    }; // class SpectralDiffer

    /**
     *
    * @class HybridDiffer
    * @brief Spectral differentiation utilities for functions defined on
    * Legendre–Gauss–Lobatto (LGL) z-nodes and uniform r-nodes.
    *
    * This class provides methods to compute spectral differentiation matrices
    * and finite difference operators
    *
    */
    class HybridDiffer {
    public:
        enum class FDOrder {
            Second,
            Fourth
        };

        HybridDiffer(size_t Nz,
                    size_t Nr,
                    Real rmin,
                    Real rmax,
                    FDOrder order = FDOrder::Fourth)
                : Nz_(Nz), Nr_(Nr), rmin_(rmin), rmax_(rmax), fd_order_(order)
        {
            if (Nr_ < 5 && fd_order_ == FDOrder::Fourth) {
                throw std::invalid_argument("FiniteDiffR: Nr must be at least 5 for 4th-order FD.");
            }
            if (Nr_ < 3) {
                throw std::invalid_argument("FiniteDiffR: Nr must be at least 3.");
            }
            if (rmax_ <= rmin_) {
                throw std::invalid_argument("FiniteDiffR: require rmax > rmin.");
            }
            // construct the z and r nodes, and precompute weights for barycentric differentiation
            z_  = SpectralDiffer::build_legendre_gauss_lobatto_nodes(Nz_);
            r_  = HybridDiffer::build_uniform_nodes(Nr_, rmin_, rmax_);
            dz_weights_ = SpectralDiffer::build_barycentric_weights_LGL(z_);
            dr_ = (rmax_ - rmin_) / Real(Nr_ - 1);
        }

        [[nodiscard]] const RVector& lgl_nodes() const noexcept { return z_; }
        [[nodiscard]] const RVector& r_nodes()   const noexcept { return r_; }
        [[nodiscard]] const RVector& z_weights() const noexcept { return dz_weights_; }

        [[nodiscard]] Real rmin() const noexcept { return rmin_; }
        [[nodiscard]] Real rmax() const noexcept { return rmax_; }
        [[nodiscard]] Real dr()   const noexcept { return dr_; }
        [[nodiscard]] size_t Nz() const noexcept { return Nz_; }
        [[nodiscard]] size_t Nr() const noexcept { return Nr_; }

        // core numerics
        void dz_barycentric_inplace(std::span<const GHPScalar<Complex>> f_span,
                                    std::span<GHPScalar<Complex>> df_dz) const;
        // wrapper
        void dz_barycentric_RSlice(const SpectralGHPVectorized::RSlice& in,
                                   SpectralGHPVectorized::RSlice& out) const;

        // core numerical kernels
        void dr_FD_inplace(std::span<const GHPScalar<Complex>> f_span,
                           std::span<GHPScalar<Complex>> df_dr) const;

        void d2r_FD_inplace(std::span<const GHPScalar<Complex>> f_span,
                            std::span<GHPScalar<Complex>> d2f_dr2) const;
        // wrappers
        void dr_FD_ZSlice(const SpectralGHPVectorized::ZSlice& in,
                          SpectralGHPVectorized::ZSlice& out) const;

        void d2r_FD_ZSlice(const SpectralGHPVectorized::ZSlice& in,
                           SpectralGHPVectorized::ZSlice& out) const;

        //static RVector build_legendre_gauss_lobatto_nodes(size_t N);
        static RVector build_uniform_nodes(size_t N, Real xmin, Real xmax);

    private:
        Real rmin_{0}, rmax_{0}, dr_{0};
        size_t Nz_{0}, Nr_{0};
        FDOrder fd_order_{FDOrder::Fourth};

        std::vector<Real> z_;
        std::vector<Real> r_;
        std::vector<Real> dz_weights_;
    };

} // spectral

#endif //GHZ_NUMERIC_SPECTRALDIFFER_HPP
