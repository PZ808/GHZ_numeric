//
// Created by Peter Zimmerman on 26.10.25.
//
#ifndef GHZ_NUMERIC_SPECTRALSOLVERS2_HPP
#define GHZ_NUMERIC_SPECTRALSOLVERS2_HPP
#pragma once

#include <vector>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "GhzTypes.hpp"
#include "GHPScalars.hpp"


namespace SpecS2boost {

    using Complex = teuk::Complex;
    using Real = teuk::Real;
    using std::vector;
    using std::pair;
    using boost::numeric::ublas::matrix;

/** \class  SpectralSolver
 *
 * \brief  Spectral solver for the 2-sphere z=(-1,1) and φ=(0,2pi)
 * constructor builds the Legendre collocation points
 * the barycenter weighs and diff matrix
 *
 */
    class SpectralSolver {
    private:
        int Nz_;       // Number of LGL nodes (θ)
        int Nphi_;     // Number of Fourier modes (φ)
        bool use_fast_double_fft_; // whether to employ double FFT

        vector<Real> z_;      // LGL nodes [-1,1]
        vector<Real> w_;      // LGL quadrature weights
        matrix<Real> D_;            // Legendre differentiation matrix

    public:
        // Constructor: build nodes, weights, D
        SpectralSolver(int Nz, int Nphi, bool use_fast_double_fft);

        // Getters
        const vector<Real>& nodes() const { return z_; }  // nodes
        const vector<Real>& weights() const { return w_; } // weights
        const matrix<Real>& Dmatrix() const { return D_; }     // diff matrix
        //const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> & Dmatrix() const { return D_; }

        // Derivatives
        // z-derivative via Gauss-Lobatto
        /**
         * \function dz
         * \brief compute the polar angular derivative using Gauss-legendre weights
         */
        void dz(const vector<vector<Complex>>& f, vector<vector<Complex>>& df_dz) const;  // polar deriv in z=cos(th)
        /**
         *  \function dphi_fft
         *   \brief φ-derivative via FFT (non-destructive)
         */
        void dphi_fft(const vector<vector<Complex>>& f,
                      vector<vector<Complex>>& df_dphi) const;

        // Barycentric weights & interpolation
        vector<Real> barycentric_weights() const;
        /**
         * \function barycentric_interp_and_derivative
         * @param f :  input function
         * @param w  : weights
         * @param z0 : point to compute the function and its derivative
         * @return pair (f(z0),f'(z0))
         * \brief interpolation routine for f and f' at z=z0
         */
        std::pair<Complex, Complex> barycentric_interp_and_derivative(
                const vector<Complex>& f, const vector<Real>& w, Real z0) const;

        // Utilities
        static std::pair<Real,Real> legendre_P_and_dP(int n, Real x);
        static vector<Real> legendre_gauss_lobatto(int N);
        static matrix<Real> legendre_diff_matrix(const vector<Real>& x);

        GHPField edth(const GHPField &f_in) const;
    };

} // namespace SpecS2


#endif //GHZ_NUMERIC_SPECTRALSOLVERS2_HPP
