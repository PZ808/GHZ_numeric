//
// Created by Peter Zimmerman on 26.10.25.
//

#ifndef GHZ_NUMERIC_SPECTRALSOLVERS2_HPP
#define GHZ_NUMERIC_SPECTRALSOLVERS2_HPP
#pragma once

#include <vector>
#include <complex>
#include "TeukTypes.hpp"
#include <Eigen/Dense>
#include <fftw3.h>

namespace SpecS2 {

    using Complex = teuk::Complex;
    using Real = teuk::Real;
    using std::vector;
    using Eigen::MatrixXd;

// ===========================================================
// Spectral solver for the 2-sphere z=(-1,1) and φ=(0,2pi)
// ===========================================================
    class SpectralSolver {
    private:
        int Nz_;       // Number of LGL nodes (θ)
        int Nphi_;     // Number of Fourier modes (φ)

        vector<Real> z_;      // LGL nodes [-1,1]
        vector<Real> w_;      // LGL quadrature weights
        MatrixXd D_;            // Legendre differentiation matrix

    public:
        // Constructor: build nodes, weights, D
        SpectralSolver(int Nz, int Nphi);

        // Getters
        const vector<Real>& nodes() const { return z_; }
        const vector<Real>& weights() const { return w_; }
        const MatrixXd& Dmatrix() const { return D_; }

        // Derivatives
        // z-derivative via Gauss-Lobatto
        void dz(const vector<vector<Complex>>& f, vector<vector<Complex>>& df_dz) const;  // polar deriv in z=cos(th)
        // φ-derivative via FFT (non-destructive)
        void dphi_fft(const vector<vector<Complex>>& f,
                      vector<vector<Complex>>& df_dphi) const;

        // Barycentric weights & interpolation
        vector<Real> barycentric_weights() const;
        std::pair<Complex, Complex> barycentric_interp_and_derivative(
                const vector<Complex>& f, const vector<Real>& w, Real z0) const;

        // Utilities
        static std::pair<Real,Real> legendre_P_and_dP(int n, Real x);
        static vector<Real> legendre_gauss_lobatto(int N);
        static MatrixXd legendre_diff_matrix(const vector<Real>& x);
    };

} // namespace SpecS2


#endif //GHZ_NUMERIC_SPECTRALSOLVERS2_HPP
