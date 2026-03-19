// Created by Peter Zimmerman on 19.03.26.
#include "ghz/orbit/KerrOrbit.hpp"
#include "ghz/geom/KerrMetric.hpp"
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>
#include <algorithm>

namespace {

    using Real = teuk::Real;
    using Complex = teuk::Complex;

    Real X_exact(Real q) {
        return std::cos(q) + Real(0.4) * std::sin(Real(2.0) * q);
    }

// Match the current sign convention in compute_Delta_and_freq_modes_SIMD
    Real Delta_exact_zero_mean(Real q) {
        return -std::sin(q) + Real(0.2) * std::cos(Real(2.0) * q);
    }

    bool approx_equal(Real a, Real b, Real tol) {
        return std::abs(a - b) < tol;
    }

    template <typename Vec>
    Real mean_of(const Vec& v) {
        if (v.empty()) return Real(0);
        Real s = std::accumulate(v.begin(), v.end(), Real(0));
        return s / Real(v.size());
    }

    Real run_nonuniform_pipeline_error(orbit::KerrBoundOrbit& orb, std::size_t N, bool verbose = true) {
        std::vector<Real> q_of_psi(N);
        std::vector<Complex> samples_in(N);

        // Build a monotone but nonuniform map q(psi) on a uniform psi-grid.
        for (std::size_t i = 0; i < N; ++i) {
            Real psi = teuk::twoPi * Real(i) / Real(N);

            // q = psi + eps sin(psi), with small eps to preserve monotonicity
            Real q = psi + Real(0.15) * std::sin(psi);
            if (q < 0) q += teuk::twoPi;
            if (q >= teuk::twoPi) q -= teuk::twoPi;
            q_of_psi[i] = q;
        }

        // monotone check
        for (std::size_t i = 1; i < N; ++i) {
            assert(q_of_psi[i] > q_of_psi[i - 1]);
        }

        // Sample X(q(psi_i))
        for (std::size_t i = 0; i < N; ++i) {
            samples_in[i] = Complex(X_exact(q_of_psi[i]), 0.0);
        }

        std::vector<Real> Delta_out;
        std::vector<Complex> modes_out;

        // Ups = 1 for the test
        orb.compute_Delta_and_freq_modes_SIMD(samples_in, q_of_psi, Real(1.0), Delta_out, modes_out);

        // Zero-mean check on reconstructed oscillatory piece
        Real mean_delta = mean_of(Delta_out);
        if (verbose) {
            std::cout << "mean(Delta_out) = " << mean_delta << "\n";
        }
        assert(std::abs(mean_delta) < 1e-10L);

        // Compare on the uniform q-grid q_j = 2pi j/N
        Real max_err = 0.0;
        for (std::size_t j = 0; j < N; ++j) {
            Real q = teuk::twoPi * Real(j) / Real(N);
            Real expected = Delta_exact_zero_mean(q);
            Real err = std::abs(Delta_out[j] - expected);
            max_err = std::max(max_err, err);
        }

        if (verbose) {
            std::cout << "N = " << N << ", nonuniform-q pipeline max error = " << max_err << "\n";
        }
        return max_err;
    }

    void test_nonuniform_q_pipeline_loose(orbit::KerrBoundOrbit& orb, std::size_t N) {
        Real max_err = run_nonuniform_pipeline_error(orb, N, true);

        // Current spline-resampling pipeline is not spectral, so keep this modest.
        assert(max_err < 5e-5L);
    }
    void test_nonuniform_q_pipeline_third_order_convergence(const KerrMetric& gKerr) {
        const Real p   = 8.0L;
        const Real e   = 0.2L;
        const Real inc = 0.4L;
        const int  chi = +1;

        const std::size_t N1 = 33;
        const std::size_t N2 = 65;
        const std::size_t N3 = 129;

        orbit::KerrBoundOrbit orb1(gKerr, p, e, inc, chi, N1, N1);
        orbit::KerrBoundOrbit orb2(gKerr, p, e, inc, chi, N2, N2);
        orbit::KerrBoundOrbit orb3(gKerr, p, e, inc, chi, N3, N3);

        orb1.initialize_orbit();
        orb2.initialize_orbit();
        orb3.initialize_orbit();

        Real e1 = run_nonuniform_pipeline_error(orb1, N1, true);
        Real e2 = run_nonuniform_pipeline_error(orb2, N2, true);
        Real e3 = run_nonuniform_pipeline_error(orb3, N3, true);

        Real h1 = Real(1.0) / Real(N1);
        Real h2 = Real(1.0) / Real(N2);
        Real h3 = Real(1.0) / Real(N3);

        Real p12 = std::log(e1 / e2) / std::log(h1 / h2);
        Real p23 = std::log(e2 / e3) / std::log(h2 / h3);

        std::cout << "Observed order p12 = " << p12 << "\n";
        std::cout << "Observed order p23 = " << p23 << "\n";

        assert(p12 > 2.8L);
        assert(p23 > 2.8L);
    }

} // namespace

int main() {
    try {
        std::cout << "Running bound orbit tests...\n";

        // Adjust this constructor line if your KerrMetric constructor differs.
        Real mass = 1.0;
        Real spin = 0.2;
        KerrParams params(mass, spin);      // set kerr params  r_+ = 1
        KerrMetric gKerr(params);

        const Real p   = 8.0L;
        const Real e   = 0.2L;
        const Real inc = 0.4L;
        const int  chi = +1;

        // Any odd values are fine for constructing the orbit object.
        const std::size_t N = 33;

        orbit::KerrBoundOrbit orb(gKerr, p, e, inc, chi, N, N);
        orb.initialize_orbit();

        test_nonuniform_q_pipeline_loose(orb, N);
        test_nonuniform_q_pipeline_third_order_convergence(gKerr);

        std::cout << "All bound orbit tests passed.\n";
        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "Test failed with exception: " << ex.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "Test failed with unknown exception.\n";
        return 1;
    }
}