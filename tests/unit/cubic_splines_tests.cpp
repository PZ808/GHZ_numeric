//
// Created by Peter Zimmerman on 19.03.26.
//
#include "ghz/orbit/Splines.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

using Real = teuk::Real;
using orbit::PeriodicSpline;

namespace {

    constexpr Real pi = Real(M_PI);
    constexpr Real twoPi = Real(2.0) * pi;

    Real f(Real x) {
        return std::sin(x) + Real(0.3) * std::cos(Real(2.0) * x);
    }

    bool approx_equal(Real a, Real b, Real tol) {
        return std::abs(a - b) < tol;
    }

    std::vector<Real> make_uniform_periodic_grid(std::size_t N) {
        std::vector<Real> x(N);
        for (std::size_t i = 0; i < N; ++i) {
            x[i] = twoPi * Real(i) / Real(N);
        }
        return x;
    }

    std::vector<Real> sample_function(const std::vector<Real>& x) {
        std::vector<Real> y(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) y[i] = f(x[i]);
        return y;
    }

    void test_knot_reproduction() {
        std::size_t N = 32;
        auto x = make_uniform_periodic_grid(N);
        auto y = sample_function(x);

        PeriodicSpline s(x, y);

        for (std::size_t i = 0; i < N; ++i) {
            assert(approx_equal(s.eval(x[i]), y[i], 1e-12L));
        }

        std::cout << "PASS: knot reproduction\n";
    }

    void test_period_wrapping() {
        std::size_t N = 32;
        auto x = make_uniform_periodic_grid(N);
        auto y = sample_function(x);

        PeriodicSpline s(x, y);

        std::vector<Real> probes = {
                Real(0.1), Real(1.3), Real(2.7), Real(5.9)
        };

        for (Real xp : probes) {
            Real a = s.eval(xp);
            Real b = s.eval(xp + twoPi);
            if (!approx_equal(a, b, 1e-10L)) {
                std::cerr << "FAIL wrapping at x=" << xp
                          << "  s(x)=" << a
                          << "  s(x+2pi)=" << b << "\n";
                assert(false);
            }
        }

        std::cout << "PASS: period wrapping\n";
    }

    void test_midpoint_accuracy() {
        std::size_t N = 64;
        auto x = make_uniform_periodic_grid(N);
        auto y = sample_function(x);

        PeriodicSpline s(x, y);

        Real max_err = 0;
        for (std::size_t i = 0; i < N; ++i) {
            Real x0 = x[i];
            Real x1 = (i + 1 < N) ? x[i + 1] : twoPi;
            Real xm = Real(0.5) * (x0 + x1);
            if (i + 1 == N) xm = Real(0.5) * (x[i] + twoPi);

            Real err = std::abs(s.eval(xm) - f(xm));
            max_err = std::max(max_err, err);
        }

        std::cout << "midpoint max err = " << max_err << "\n";
        assert(max_err < 1e-4L); // conservative
        std::cout << "PASS: midpoint accuracy\n";
    }

    void test_near_endpoint_consistency() {
        std::size_t N = 64;
        auto x = make_uniform_periodic_grid(N);
        auto y = sample_function(x);

        PeriodicSpline s(x, y);

        Real eps = 1e-6L;
        Real a = s.eval(eps);
        Real b = s.eval(twoPi + eps);
        Real c = s.eval(-eps);
        Real d = s.eval(twoPi - eps);

        assert(approx_equal(a, b, 1e-10L));
        assert(approx_equal(c, d, 1e-10L));

        std::cout << "PASS: near-endpoint consistency\n";
    }

    void test_constant_function() {
        std::size_t N = 16;
        auto x = make_uniform_periodic_grid(N);
        std::vector<Real> y(N, Real(2.5));

        PeriodicSpline s(x, y);

        for (Real xp : {Real(-1.0), Real(0.0), Real(1.0), Real(7.0), Real(10.0)}) {
            assert(approx_equal(s.eval(xp), Real(2.5), 1e-14L));
        }

        std::cout << "PASS: constant function\n";
    }

} // namespace

int main() {
    test_constant_function();
    test_knot_reproduction();
    test_period_wrapping();
    test_midpoint_accuracy();
    test_near_endpoint_consistency();

    std::cout << "All PeriodicSpline tests passed.\n";
    return 0;
}
