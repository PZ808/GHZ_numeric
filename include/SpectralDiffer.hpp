//
// Created by Peter Zimmerman on 28.11.25.
//

#ifndef GHZ_NUMERIC_SPECTRALDIFFER_HPP
#define GHZ_NUMERIC_SPECTRALDIFFER_HPP

#include "GhzTypes.hpp"
#include "HeldScalars.hpp"
#include "VectorsGHZ.hpp"  // Your fixed-size vector class
#include "SpectralGHPField.hpp"
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

    /**
     * @class SpectralDiffer
     * @brief Spectral differentiation utilities for functions defined on
     * Legendre–Gauss–Lobatto (LGL) and Chebyshev–Lobatto (CL) nodes.
     *
     * This class provides methods to compute spectral differentiation matrices
     * and perform differentiation operations on functions sampled at LGL and CL nodes.
     * It includes routines for constructing compu
     *  edth/edthbar operators, as well as thorn
     *
     */
    class SpectralDiffer {

    public:
        // Constructor: build nodes, weights, differentiation matrix
        SpectralDiffer(size_t Nz, size_t Nr) : Nz_(Nz), Nr_(Nr) {
            z_ = build_legendre_gauss_lobatto_nodes(Nz_); // generate nodes z[i]
            r_ = build_chebyshev_lobatto_nodes(Nr_); // generate nodes z[i]
            wz_ = build_barycentric_weights_from_nodes(z_);
            wr_ = build_barycentric_weights_from_nodes(r_);
            w_ = build_barycentric_weights_z();
            Dz_ = build_legendre_diff_matrix(z_);
            Dr_ = build_chebyshev_diff_matrix(r_);
        }

        // Getters
        [[nodiscard]] const RVector &lgl_nodes() const { return z_; }
        [[nodiscard]] const RVector &cl_nodes() const { return r_; }
        [[nodiscard]] const RVector &z_weights() const { return wz_; }
        [[nodiscard]] const RVector &r_weights() const { return wr_; }

        [[nodiscard]] const matrix<Real> &Dz_matrix() const { return Dz_; }
        [[nodiscard]] const matrix<Real> &Dr_matrix() const { return Dr_; }

        // Barycentric utilities (declare; implement elsewhere or below)
        [[nodiscard]] RVector build_barycentric_weights_z() const;
        RVector build_barycentric_weights_from_nodes(const RVector& nodes) const;

        [[nodiscard]] static std::pair<Complex, Complex> barycentric_interp_and_derivative(const RVector& nodes,
                                                                                           const CVector& f,
                                                                                           const RVector& w,
                                                                                           Real x0);

        static std::pair<Real, Real> legendre_P_and_dP(unsigned int n, const Real &x);

        static std::vector<Real> build_legendre_gauss_lobatto_nodes(int N);

        static std::vector<Real> build_chebyshev_lobatto_nodes(int N);

        static matrix<Real> build_legendre_diff_matrix(const RVector& x);
        static matrix<Real> build_chebyshev_diff_matrix(const RVector& x);


        // dz on a single r=const slice (in-place into df_dz)
        void dz_Dmatrix_RSlice(const GHPSpectralField::RSlice &f_RSlice,
                               GHPSpectralField::RSlice &dfdz_RSlice_into) const;

        void dr_Dmatrix_ZSlice(const  spectral::GHPSpectralField::ZSlice &f_ZSlice,
                               GHPSpectralField::ZSlice &df_dr) const;

        void dz_Dmatrix(const std::span<GHPScalar<Complex>> &f_span_rconst,
                        std::span<GHPScalar<Complex>> &df_dz) const;

        void dr_Dmatrix(const std::span<GHPScalar<Complex>>& f_span_zconst,
                        std::span<GHPScalar<Complex>>& df_dr) const;


        GHPSpectralField::RSlice dz_barycentric_RSlice(const GHPSpectralField::RSlice& f_RSlice,
                                                       GHPSpectralField::RSlice& dfdz_RSlice_into,
                                                       const std::vector<Real>& weights) const;

        void dz_barycentric_inplace(const std::span<GHPScalar<Complex>> &f_span_const,
                                    std::span<GHPScalar<Complex>> &df,
                                    const std::vector<Real> &weights) const;

        GHPSpectralField::ZSlice dr_barycentric_ZSlice(const GHPSpectralField::ZSlice &f,
                                                      GHPSpectralField::ZSlice &df_dr,
                                                      const std::vector<Real> &w ) const;

        void dr_barycentric_inplace(const std::span<GHPScalar<Complex>> &f_span_const,
                                    std::span<GHPScalar<Complex>> &df,
                                    const std::vector<Real> &weights) const;

        void dphi_single_m_slice(const std::span<GHPScalar<Complex>> &f_slice,
                                 std::span<GHPScalar<Complex>> &df_dphi_slice,
                                 const Real &m) const;

        void dt_single_omega_slice(const std::span<GHPScalar<Complex>> &f_slice,
                                   std::span<GHPScalar<Complex>> &df_dt_slice,
                                   const Real &omega) const;

        void edthH_bary_inplace(const GHPSpectralField::RSlice &f,
                                                GHPSpectralField::RSlice &out, Real const& a) const;
        void edthH_dmat_RSlice_inplace(const GHPSpectralField::RSlice &f,
                                       GHPSpectralField::RSlice &out, const Real& a) const;

        void edthH_inplace_RSliceV(
                const SpectralGHPVectorized::RSlice &in_RSlice,
                SpectralGHPVectorized::RSlice &out_RSlice, const Real& a) const;

        void edthBarH_inplace_RSliceV(
                const SpectralGHPVectorized::RSlice &in_RSlice,
                SpectralGHPVectorized::RSlice &out_RSlice,
                const Real& a) const;

    private:
        size_t Nz_, Nr_;                  // Number of nodes
        std::vector<Real> z_;            // LGL nodes [-1,1]
        std::vector<Real> r_;            // CL nodes [-1,1]
        std::vector<Real> w_, wz_, wr_;  // barycentric weights
        matrix<Real> Dz_;                // Legendre differentiation matrix
        matrix<Real> Dr_;                // Chebyshev differentiation matrix
    };

} // spectral
#endif //GHZ_NUMERIC_SPECTRALDIFFER_HPP
