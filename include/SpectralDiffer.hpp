//
// Created by Peter Zimmerman on 28.11.25.
//

#ifndef GHZ_NUMERIC_SPECTRALDIFFER_HPP
#define GHZ_NUMERIC_SPECTRALDIFFER_HPP

#include "GhzTypes.hpp"
#include "HeldScalars.hpp"
#include "VectorsGHZ.hpp"  // Your fixed-size vector class
#include "SpectralGHPField.hpp"

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

namespace  spectral {

    using Complex = teuk::Complex;
    using Real = teuk::Real;
    using boost::numeric::ublas::matrix;
// -----------------------------------------------------------------------------
// SpectralDiffer: operates on Z-slices for a given block (fixed w,m)
// -----------------------------------------------------------------------------
    class SpectralDiffer {
    private:
        int Nz_, Nr_;                    // Number of nodes
        bool use_fast_double_fft_;      // Whether to use double FFT

        std::vector <Real> z_;           // LGL nodes [-1,1]
        std::vector <Real> r_;           // CL nodes [-1,1]
        std::vector <Real> w_;           // LGL quadrature weights
        matrix <Real> D_;                // Legendre differentiation matrix

    public:
        // Constructor: build nodes, weights, differentiation matrix
        SpectralDiffer(int Nz, int Nr, bool use_fast_double_fft = true)
                : Nz_(Nz), Nr_(Nr), use_fast_double_fft_(use_fast_double_fft) {
            z_ = build_legendre_gauss_lobatto_nodes(Nz_); // generate nodes z[i]
            r_ = build_chebyshev_lobatto_nodes(Nr_); // generate nodes z[i]
            w_ = build_barycentric_weights();
            D_ = build_legendre_diff_matrix(z_);
        }

        // Getters
        const std::vector <Real> &lgl_nodes() const { return z_; }
        const std::vector <Real> &cl_nodes() const { return r_; }

        const std::vector <Real> &weights() const { return w_; }

        const matrix <Real> &Dmatrix() const { return D_; }

        // Barycentric utilities (declare; implement elsewhere or below)
        std::vector <Real> build_barycentric_weights() const;

        std::pair <Complex, Complex> barycentric_interp_and_derivative(
                const std::vector <Complex> &f, const std::vector <Real> &w, Real z0) const;

        static std::pair <Real, Real> legendre_P_and_dP(int n, Real x);

        static std::vector <Real> build_legendre_gauss_lobatto_nodes(int N);
        static std::vector <Real> build_chebyshev_lobatto_nodes(int N);

        static matrix <Real> build_legendre_diff_matrix(const std::vector <Real> &x);

        // dz on a single ZSlice (in-place into df_dz)
        void dz_Dmatrix(const GHPSpectralField::ZSlice &f, GHPSpectralField::ZSlice &df_dz) const;

        GHPSpectralField::ZSlice dz_barycentric(const GHPSpectralField::ZSlice &f,
                                                GHPSpectralField::ZSlice &df_dz,
                                                const std::vector <Real> &w) const;

        // d/dphi for a single slice
        void dphi_fft_single_m(const GHPSpectralField::ZSlice &f_slice,
                               GHPSpectralField::ZSlice &df_dphi_slice) const;
        void dphi_fft_single_w(const GHPSpectralField::ZSlice &f_slice,
                               GHPSpectralField::ZSlice &df_dphi_slice) const;

        // edth operator on a single slice -> returns a new ZSlice-owned buffer (or you can provide out)
        GHPSpectralField::ZSlice edth(const GHPSpectralField::ZSlice &f) const;




        // ------------------------------
        // Differentiation along z
        // ------------------------------

        /**
         * @brief Compute dz using the differentiation matrix for a single row.
         *
         * @param f_slice Input slice (row of GHPFieldVectorized)
         * @param df_dz  Output slice (preallocated)
         */
        void dz_Dmatrix_slice(const std::span<GHPScalar<Complex>>& f_slice,
                              std::span<GHPScalar<Complex>>& df_dz) const;

        /**
         * @brief Compute dz using barycentric differentiation.
         *
         * @param f_slice Input slice
         * @param df_dz  Output slice (preallocated)
         * @param w       Barycentric weights
         */
        void dz_barycentric_slice(const std::span<GHPScalar<Complex>>& f_slice,
                            std::span<GHPScalar<Complex>>& df_dz,
                            const std::vector<Real>& w) const;

        // ------------------------------
        // Derivatives along phi
        // ------------------------------

        /**
         * @brief Compute d/dphi using FFT for a single m-mode.
         *
         * @param f_slice Input slice
         * @param df_dphi_slice Output slice (preallocated)
         */
        void dphi_fft_single_m_slice(const std::span<GHPScalar<Complex>>& f_slice,
                               std::span<GHPScalar<Complex>>& df_dphi_slice,
                               const int& m) const;

        void dphi_fft_single_w_slice(const std::span<GHPScalar<Complex>>& f_slice,
                               std::span<GHPScalar<Complex>>& df_dphi_slice,
                               const int& m) const;

        /**
         * @brief Apply edth operator to a single slice.
         *
         * Returns a new temporary vector (or you can provide preallocated slice).
         */
        std::vector<GHPScalar<Complex>> edth_slice(const std::span<GHPScalar<Complex>>& f_slice, int m) const;
    };

}; // spectral

#endif //GHZ_NUMERIC_SPECTRALDIFFER_HPP
