//
// Created by Peter Zimmerman on 07.03.26.
//
// spectral_unit_tests.cpp
//
// Basic unit tests for spectral::SpectralDiffer and spectral::HybridDiffer.
// Focuses on:
//   - node generation sanity
//   - barycentric interpolation/derivative
//   - z-derivatives on LGL nodes (D-matrix + barycentric)
//   - r-derivatives on Chebyshev-Lobatto nodes (D-matrix + barycentric)
//   - r-derivatives on uniform nodes (FD first + second derivatives)
//
// Compile example:
//   c++ -O2 -std=c++20 -fopenmp spectral_unit_tests.cpp SpectralDiffer.cpp -I/path/to/include -o unit_tests
//
//

#include "ghz/spectral/SpectralDiffer.hpp"

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

namespace {

    using spectral::SpectralDiffer;
    using spectral::HybridDiffer;
    using spectral::RVector;
    using teuk::Real;
    using teuk::Complex;

    int p_test = 1;
    int q_test = -3;

    template <typename T>
    using GHP = GHPScalar<T>;

    struct TestFailure : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    inline GHP<Complex> make_scalar(const Complex& v, int p = 0, int q = 0) {
        return GHP<Complex>(v, p, q);
    }

    inline Real abs_err(const Complex& a, const Complex& b) {
        return std::abs(a - b);
    }

    inline void require_true(bool cond, const std::string& msg) {
        if (!cond) throw TestFailure(msg);
    }

    inline void require_near(const Complex& got, const Complex& expected,
                             Real tol, const std::string& msg) {
        const Real err = abs_err(got, expected);
        if (err > tol) {
            std::ostringstream oss;
            oss << msg << "\n"
                << "  got      = " << got << "\n"
                << "  expected = " << expected << "\n"
                << "  |err|    = " << err << "\n"
                << "  tol      = " << tol;
            throw TestFailure(oss.str());
        }
    }

    inline void require_near_real(Real got, Real expected,
                                  Real tol, const std::string& msg) {
        const Real err = std::abs(got - expected);
        if (err > tol) {
            std::ostringstream oss;
            oss << msg << "\n"
                << "  got      = " << got << "\n"
                << "  expected = " << expected << "\n"
                << "  |err|    = " << err << "\n"
                << "  tol      = " << tol;
            throw TestFailure(oss.str());
        }
    }

    template <typename F>
    std::vector<GHP<Complex>> sample_on_nodes(const RVector& x, F f, int p = 0, int q = 0) {
        std::vector<GHP<Complex>> out(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            out[i] = make_scalar(f(x[i]), p, q);
        }
        return out;
    }

    template <typename DF>
    void check_derivative_vector(const RVector& x,
                                 const std::vector<GHP<Complex>>& got,
                                 DF df_exact,
                                 Real tol,
                                 const std::string& label) {
        require_true(got.size() == x.size(), label + ": size mismatch");
        for (size_t i = 0; i < x.size(); ++i) {
            const Complex expected = df_exact(x[i]);
            require_near(got[i].value(), expected, tol,
                         label + " at node i=" + std::to_string(i));
        }
    }

    void test_legendre_endpoints() {
        auto [P0, dP0] = SpectralDiffer::legendre_P_and_dP(0u, Real(1));
        require_near_real(P0, Real(1), Real(0), "P0(1)");
        require_near_real(dP0, Real(0), Real(0), "P0'(1)");

        auto [P1m, dP1m] = SpectralDiffer::legendre_P_and_dP(1u, Real(-1));
        require_near_real(P1m, Real(-1), Real(0), "P1(-1)");
        require_near_real(dP1m, Real(1), Real(0), "P1'(-1)");

        // Check endpoint derivative formula P'_n(1)=n(n+1)/2
        for (unsigned int n = 2; n <= 8; ++n) {
            auto [Pp, dPp] = SpectralDiffer::legendre_P_and_dP(n, Real(1));
            auto [Pm, dPm] = SpectralDiffer::legendre_P_and_dP(n, Real(-1));

            require_near_real(Pp, Real(1), Real(0), "Pn(1)");
            require_near_real(dPp, Real(n * (n + 1) / 2.0), Real(0), "Pn'(1)");

            const Real expectedPm = (n % 2 == 0) ? Real(1) : Real(-1);
            const Real expecteddPm =
                    (n % 2 == 0) ? -Real(n * (n + 1) / 2.0) : Real(n * (n + 1) / 2.0);

            require_near_real(Pm, expectedPm, Real(0), "Pn(-1)");
            require_near_real(dPm, expecteddPm, Real(0), "Pn'(-1)");
        }
    }

    void test_node_builders() {
        constexpr size_t Nz = 12;
        constexpr size_t Nr = 15;

        SpectralDiffer diff(Nz, Nr);

        const auto& z = diff.lgl_nodes();
        const auto& r = diff.cl_nodes();

        require_true(z.size() == Nz, "LGL node count");
        require_true(r.size() == Nr, "Cheb-Lobatto node count");

        require_near_real(z.front(), Real(-1), Real(1e-14), "LGL left endpoint");
        require_near_real(z.back(),  Real(1),  Real(1e-14), "LGL right endpoint");

        require_near_real(r.front(), Real(-1), Real(1e-14), "Cheb left endpoint");
        require_near_real(r.back(),  Real(1),  Real(1e-14), "Cheb right endpoint");

        for (size_t i = 1; i < z.size(); ++i) {
            require_true(z[i] > z[i-1], "LGL nodes strictly increasing");
        }
        for (size_t i = 1; i < r.size(); ++i) {
            require_true(r[i] > r[i-1], "Cheb nodes strictly increasing");
        }

        HybridDiffer hdiff(10, 21, Real(-3), Real(7), HybridDiffer::FDOrder::Fourth);
        const auto& ru = hdiff.r_nodes();
        require_true(ru.size() == 21, "uniform r count");
        require_near_real(ru.front(), Real(-3), Real(1e-14), "uniform left endpoint");
        require_near_real(ru.back(),  Real(7),  Real(1e-14), "uniform right endpoint");

        const Real h = (Real(7) - Real(-3)) / Real(20);
        for (size_t i = 1; i < ru.size(); ++i) {
            require_near_real(ru[i] - ru[i-1], h, Real(1e-14), "uniform spacing");
        }
    }
    void test_barycentric_interp_and_derivative() {
        constexpr size_t Nz = 14;
        SpectralDiffer diff(Nz, 10);

        const auto& x = diff.lgl_nodes();
        const auto  w = SpectralDiffer::build_barycentric_weights_LGL(x);

        auto f = [](Real z) -> Complex {
            return Complex(z*z*z - Real(2)*z + Real(1),
                           Real(2)*z*z + Real(3)*z - Real(4));
        };
        auto df = [](Real z) -> Complex {
            return Complex(Real(3)*z*z - Real(2),
                           Real(4)*z + Real(3));
        };

        std::vector<Complex> fv(x.size());
        for (size_t i = 0; i < x.size(); ++i) fv[i] = f(x[i]);

        std::cout << "LGL nodes and weights:\n";
        for (size_t i = 0; i < x.size(); ++i) {
            std::cout << i << "  x=" << x[i] << "  w=" << w[i] << "\n";
        }

        // First: test exactly at every node
        for (size_t i = 0; i < x.size(); ++i) {
            auto [fi, dfi] = SpectralDiffer::barycentric_interp_and_derivative(x, fv, w, x[i]);
            require_near(fi, f(x[i]), Real(1e-14),
                         "interp reproduces nodal value at i=" + std::to_string(i));
        }

        // Then: one explicit node test including derivative
        const size_t k = x.size() / 2;
        auto [fk, dfk] = SpectralDiffer::barycentric_interp_and_derivative(x, fv, w, x[k]);
        require_near(fk,  f(x[k]),  Real(1e-14), "barycentric interp at node");
        require_near(dfk, df(x[k]), Real(1e-11), "barycentric derivative at node");

        // Finally: off-grid point
        const Real z0 = Real(0.137);
        auto [f0, df0] = SpectralDiffer::barycentric_interp_and_derivative(x, fv, w, z0);

        require_near(f0,  f(z0),  Real(1e-12), "barycentric interp at off-grid point");
        require_near(df0, df(z0), Real(1e-11), "barycentric derivative at off-grid point");
    }

    void test_dz_spectral_kernels() {
        constexpr size_t Nz = 18;
        constexpr size_t Nr = 12;
        SpectralDiffer diff(Nz, Nr);

        const auto& z = diff.lgl_nodes();

        auto f = [](Real x) -> Complex {
            return Complex(x*x*x*x - Real(2)*x*x + Real(3)*x - Real(1),
                           Real(2)*x*x*x - Real(5)*x + Real(7));
        };
        auto df = [](Real x) -> Complex {
            return Complex(Real(4)*x*x*x - Real(4)*x + Real(3),
                           Real(6)*x*x - Real(5));
        };

        auto in  = sample_on_nodes(z, f, p_test, q_test);
        std::vector<GHP<Complex>> outD(Nz, make_scalar(Complex(0,0), p_test, q_test));
        std::vector<GHP<Complex>> outB(Nz, make_scalar(Complex(0,0), p_test, q_test));

        diff.dz_Dmatrix(std::span<const GHP<Complex>>(in.data(), in.size()),
                        std::span<GHP<Complex>>(outD.data(), outD.size()));

        diff.dz_barycentric_inplace(std::span<const GHP<Complex>>(in.data(), in.size()),
                                    std::span<GHP<Complex>>(outB.data(), outB.size()));

        check_derivative_vector(z, outD, df, Real(5e-12), "dz_Dmatrix polynomial exactness");
        check_derivative_vector(z, outB, df, Real(5e-12), "dz_barycentric polynomial exactness");
    }

    void test_dr_spectral_kernels() {
        constexpr size_t Nz = 10;
        constexpr size_t Nr = 20;
        SpectralDiffer diff(Nz, Nr);

        const auto& r = diff.cl_nodes();

        auto f = [](Real x) -> Complex {
            return Complex(Real(2)*x*x*x - Real(3)*x*x + Real(4)*x - Real(5),
                           -x*x*x*x + Real(2)*x*x - Real(1));
        };
        auto df = [](Real x) -> Complex {
            return Complex(Real(6)*x*x - Real(6)*x + Real(4),
                           -Real(4)*x*x*x + Real(4)*x);
        };

        auto in  = sample_on_nodes(r, f, p_test, q_test);
        std::vector<GHP<Complex>> outD(Nr, make_scalar(Complex(0,0), p_test, q_test));
        std::vector<GHP<Complex>> outB(Nr, make_scalar(Complex(0,0), p_test, q_test));

        diff.dr_Dmatrix(std::span<const GHP<Complex>>(in.data(), in.size()),
                        std::span<GHP<Complex>>(outD.data(), outD.size()));

        diff.dr_barycentric_inplace(std::span<const GHP<Complex>>(in.data(), in.size()),
                                    std::span<GHP<Complex>>(outB.data(), outB.size()));

        check_derivative_vector(r, outD, df, Real(5e-11), "dr_Dmatrix polynomial exactness");
        check_derivative_vector(r, outB, df, Real(5e-11), "dr_barycentric polynomial exactness");
    }

    void test_hybrid_fd_first_derivative_second_order() {
        constexpr size_t Nz = 10;
        constexpr size_t Nr = 17;
        HybridDiffer hdiff(Nz, Nr, Real(-2), Real(3), HybridDiffer::FDOrder::Second);
        const auto& r = hdiff.r_nodes();

        // Linear function should be differentiated exactly by both interior and one-sided stencils.
        auto f = [](Real x) -> Complex {
            return Complex(Real(3)*x - Real(2), -Real(4)*x + Real(1));
        };
        auto df = [](Real) -> Complex {
            return Complex(Real(3), -Real(4));
        };

        auto in = sample_on_nodes(r, f, p_test, q_test);
        std::vector<GHP<Complex>> out(Nr, make_scalar(Complex(0,0), p_test, q_test));

        hdiff.dr_FD_inplace(std::span<const GHP<Complex>>(in.data(), in.size()),
                            std::span<GHP<Complex>>(out.data(), out.size()));

        check_derivative_vector(r, out, df, Real(1e-13), "Hybrid FD(2) first derivative linear exactness");
    }

    void test_hybrid_fd_second_derivative_second_order() {
        constexpr size_t Nz = 10;
        constexpr size_t Nr = 17;
        HybridDiffer hdiff(Nz, Nr, Real(-1), Real(2), HybridDiffer::FDOrder::Second);
        const auto& r = hdiff.r_nodes();

        // Quadratic should give exact second derivative with these stencils.
        auto f = [](Real x) -> Complex {
            return Complex(Real(2)*x*x - Real(3)*x + Real(1),
                           -Real(5)*x*x + Real(4)*x - Real(7));
        };
        auto d2f = [](Real) -> Complex {
            return Complex(Real(4), -Real(10));
        };

        auto in = sample_on_nodes(r, f, p_test, q_test);
        std::vector<GHP<Complex>> out(Nr, make_scalar(Complex(0,0), p_test, q_test));

        hdiff.d2r_FD_inplace(std::span<const GHP<Complex>>(in.data(), in.size()),
                             std::span<GHP<Complex>>(out.data(), out.size()));

        check_derivative_vector(r, out, d2f, Real(1e-12), "Hybrid FD(2) second derivative quadratic exactness");
    }

    void test_hybrid_fd_first_derivative_fourth_order() {
        constexpr size_t Nz = 10;
        constexpr size_t Nr = 21;
        HybridDiffer hdiff(Nz, Nr, Real(-2), Real(2), HybridDiffer::FDOrder::Fourth);
        const auto& r = hdiff.r_nodes();

        // Cubic should be exact for fourth-order first-derivative stencils.
        auto f = [](Real x) -> Complex {
            return Complex(x*x*x - Real(2)*x*x + Real(5)*x - Real(3),
                           -Real(2)*x*x*x + x*x - Real(7)*x + Real(1));
        };
        auto df = [](Real x) -> Complex {
            return Complex(Real(3)*x*x - Real(4)*x + Real(5),
                           -Real(6)*x*x + Real(2)*x - Real(7));
        };

        auto in = sample_on_nodes(r, f, p_test,q_test);
        std::vector<GHP<Complex>> out(Nr, make_scalar(Complex(0,0), p_test, q_test));

        hdiff.dr_FD_inplace(std::span<const GHP<Complex>>(in.data(), in.size()),
                            std::span<GHP<Complex>>(out.data(), out.size()));

        check_derivative_vector(r, out, df, Real(1e-11), "Hybrid FD(4) first derivative cubic exactness");
    }

    void test_hybrid_fd_second_derivative_fourth_order() {
        constexpr size_t Nz = 10;
        constexpr size_t Nr = 21;
        HybridDiffer hdiff(Nz, Nr, Real(-3), Real(1), HybridDiffer::FDOrder::Fourth);
        const auto& r = hdiff.r_nodes();

        // Quartic should be exact for fourth-order second-derivative stencils.
        auto f = [](Real x) -> Complex {
            return Complex(x*x*x*x - Real(2)*x*x*x + Real(3)*x*x - Real(5)*x + Real(2),
                           -Real(2)*x*x*x*x + x*x*x - Real(4)*x*x + Real(6));
        };
        auto d2f = [](Real x) -> Complex {
            return Complex(Real(12)*x*x - Real(12)*x + Real(6),
                           -Real(24)*x*x + Real(6)*x - Real(8));
        };

        auto in = sample_on_nodes(r, f, p_test, q_test);
        std::vector<GHP<Complex>> out(Nr, make_scalar(Complex(0,0), p_test, q_test));

        hdiff.d2r_FD_inplace(std::span<const GHP<Complex>>(in.data(), in.size()),
                             std::span<GHP<Complex>>(out.data(), out.size()));

        check_derivative_vector(r, out, d2f, Real(1e-10), "Hybrid FD(4) second derivative quartic exactness");
    }

    void test_hybrid_dz_kernel() {
        constexpr size_t Nz = 16;
        constexpr size_t Nr = 11;
        HybridDiffer hdiff(Nz, Nr, Real(-1), Real(1), HybridDiffer::FDOrder::Fourth);
        const auto& z = hdiff.lgl_nodes();

        auto f = [](Real x) -> Complex {
            return Complex(-x*x*x + Real(2)*x*x + Real(3),
                           Real(5)*x*x*x*x - x + Real(2));
        };
        auto df = [](Real x) -> Complex {
            return Complex(-Real(3)*x*x + Real(4)*x,
                           Real(20)*x*x*x - Real(1));
        };

        auto in = sample_on_nodes(z, f, p_test, q_test);
        std::vector<GHP<Complex>> out(Nz, make_scalar(Complex(0,0), p_test, q_test));

        hdiff.dz_barycentric_inplace(std::span<const GHP<Complex>>(in.data(), in.size()),
                                     std::span<GHP<Complex>>(out.data(), out.size()));

        check_derivative_vector(z, out, df, Real(5e-12), "Hybrid dz barycentric exactness");
    }

    using TestFn = void(*)();

    struct TestCase {
        const char* name;
        TestFn fn;
    };

} // namespace

int main() {
    const std::vector<TestCase> tests = {
            {"legendre_endpoints", test_legendre_endpoints},
            {"node_builders", test_node_builders},
            {"barycentric_interp_and_derivative", test_barycentric_interp_and_derivative},
            {"dz_spectral_kernels", test_dz_spectral_kernels},
            {"dr_spectral_kernels", test_dr_spectral_kernels},
            {"hybrid_fd_first_derivative_second_order", test_hybrid_fd_first_derivative_second_order},
            {"hybrid_fd_second_derivative_second_order", test_hybrid_fd_second_derivative_second_order},
            {"hybrid_fd_first_derivative_fourth_order", test_hybrid_fd_first_derivative_fourth_order},
            {"hybrid_fd_second_derivative_fourth_order", test_hybrid_fd_second_derivative_fourth_order},
            {"hybrid_dz_kernel", test_hybrid_dz_kernel},
    };

    int passed = 0;
    for (const auto& t : tests) {
        try {
            t.fn();
            std::cout << "[PASS] " << t.name << "\n";
            ++passed;
        } catch (const std::exception& e) {
            std::cerr << "[FAIL] " << t.name << "\n" << e.what() << "\n";
            return EXIT_FAILURE;
        } catch (...) {
            std::cerr << "[FAIL] " << t.name << "\nUnknown exception\n";
            return EXIT_FAILURE;
        }
    }

    std::cout << "\nAll " << passed << " tests passed.\n";
    return EXIT_SUCCESS;
}