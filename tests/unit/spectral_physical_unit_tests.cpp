//
// Created by Peter Zimmerman on 18.03.26.
//
//
// spectral_physical_unit_tests.cpp
//
// Focused tests for the physical radial layer built from
//   - spectral::SpectralDiffer   (reference x-grid)
//   - ghz::numeric::AffineMap1D  (x <-> r)
//   - ghz::numeric::PhysicalChebRadialOps
//
// What is tested:
//   1. Physical Dr and Drr acting on simple functions of physical r
//   2. Off-grid interpolation and physical derivative evaluation
//   3. Manufactured one-domain Schwarzschild collocation solve
//
// Compile example:
//   c++ -O2 -std=c++20 -fopenmp spectral_physical_unit_tests.cpp \
//      SpectralDiffer.cpp -I/path/to/include -o spectral_physical_unit_tests
//

#include "ghz/core/GhzTypes.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/spectral/SpectralCoordinateMaps.hpp"
#include "ghz/spectral/PhysicalChebRadialOps.hpp"
#include "ghz/geom/Domain.hpp"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

    using teuk::Real;
    using teuk::Complex;

    using spectral::SpectralDiffer;
    using spectral::matrix;
    using spectral::RVector;
    using spectral::CVector;

    namespace ublas = boost::numeric::ublas;

    // -------------------------------------------------------------------------
    // Simple testing helpers
    // -------------------------------------------------------------------------

    [[noreturn]] void fail(const std::string& msg) {
        throw std::runtime_error(msg);
    }

    void require_true(bool cond, const std::string& msg) {
        if (!cond) fail(msg);
    }

    void require_close_real(Real got, Real expected, Real tol, const std::string& label) {
        const Real err = std::abs(got - expected);
        if (err > tol) {
            std::ostringstream os;
            os << label << " failed: got=" << std::setprecision(18) << got
               << " expected=" << expected
               << " abs_err=" << err
               << " tol=" << tol;
            fail(os.str());
        }
    }

    void require_close_complex(Complex got, Complex expected, Real tol, const std::string& label) {
        const Real err = std::abs(got - expected);
        if (err > tol) {
            std::ostringstream os;
            os << label << " failed: got=" << std::setprecision(18) << got
               << " expected=" << expected
               << " abs_err=" << err
               << " tol=" << tol;
            fail(os.str());
        }
    }

    template <typename F, typename G>
    Real max_abs_error_on_grid(const std::vector<Real>& x, F f, G g) {
        Real err = Real(0);
        for (const auto& xi : x) {
            err = std::max(err, std::abs(f(xi) - g(xi)));
        }
        return err;
    }

    // -------------------------------------------------------------------------
    // Matrix-vector helpers
    // -------------------------------------------------------------------------

    std::vector<Real> apply_real_matrix(
            const matrix<Real>& A,
            const std::vector<Real>& x)
    {
        const std::size_t N = x.size();
        require_true(A.size1() == N && A.size2() == N,
                     "apply_real_matrix: size mismatch.");

        std::vector<Real> y(N, Real(0));
        for (std::size_t i = 0; i < N; ++i) {
            Real sum = Real(0);
            for (std::size_t j = 0; j < N; ++j) {
                sum += A(i, j) * x[j];
            }
            y[i] = sum;
        }
        return y;
    }

    std::vector<Complex> solve_real_matrix_complex_rhs(
            const matrix<Real>& A,
            const std::vector<Complex>& b)
    {
        const std::size_t N = b.size();
        require_true(A.size1() == N && A.size2() == N,
                     "solve_real_matrix_complex_rhs: size mismatch.");

        matrix<Real> LU(A);
        ublas::permutation_matrix<std::size_t> pm(N);

        const int status = ublas::lu_factorize(LU, pm);
        if (status != 0) {
            fail("solve_real_matrix_complex_rhs: LU factorization failed.");
        }

        ublas::vector<Real> bre(N), bim(N);
        for (std::size_t i = 0; i < N; ++i) {
            bre(i) = std::real(b[i]);
            bim(i) = std::imag(b[i]);
        }

        ublas::lu_substitute(LU, pm, bre);
        ublas::lu_substitute(LU, pm, bim);

        std::vector<Complex> x(N);
        for (std::size_t i = 0; i < N; ++i) {
            x[i] = Complex(bre(i), bim(i));
        }
        return x;
    }

    // -------------------------------------------------------------------------
    // Test 1: physical Dr and Drr
    // -------------------------------------------------------------------------

    void test_physical_derivatives()
    {
        std::cout << "[1/3] testing physical Dr and Drr ...\n";

        constexpr std::size_t Nz = 9;
        constexpr std::size_t Nr = 33;

        SpectralDiffer differ(Nz, Nr);

        // Adjust this type/name if your Domain namespace differs.
        ghz::numeric::Domain domain{2.0, 6.0};
        ghz::numeric::PhysicalChebRadialOps rops(differ, domain);

        const auto& r   = rops.r();
        const auto& Dr  = rops.Dr_matrix();
        const auto& Drr = rops.Drr_matrix();

        const std::size_t N = r.size();

        std::vector<Real> f_const(N), f_r(N), f_r2(N);
        for (std::size_t i = 0; i < N; ++i) {
            f_const[i] = Real(1);
            f_r[i]     = r[i];
            f_r2[i]    = r[i] * r[i];
        }

        const auto Dr_const = apply_real_matrix(Dr,  f_const);
        const auto Dr_r     = apply_real_matrix(Dr,  f_r);
        const auto Drr_r2   = apply_real_matrix(Drr, f_r2);

        Real err_const = Real(0);
        Real err_r     = Real(0);
        Real err_r2    = Real(0);

        for (std::size_t i = 0; i < N; ++i) {
            err_const = std::max(err_const, std::abs(Dr_const[i] - Real(0)));
            err_r     = std::max(err_r,     std::abs(Dr_r[i]     - Real(1)));
            err_r2    = std::max(err_r2,    std::abs(Drr_r2[i]   - Real(2)));
        }

        std::cout << "    max |Dr[1]|      = " << std::setprecision(18) << err_const << "\n";
        std::cout << "    max |Dr[r]-1|    = " << std::setprecision(18) << err_r     << "\n";
        std::cout << "    max |Drr[r^2]-2| = " << std::setprecision(18) << err_r2    << "\n";

        const Real eps = std::numeric_limits<Real>::epsilon();
        const Real tol_drr = Real(50) * eps * std::pow(Real(N), 4);

        require_true(err_const < Real(1e-11), "physical Dr[1]=0 test failed");
        require_true(err_r     < Real(1e-11), "physical Dr[r]=1 test failed");
        require_true(err_r2    < tol_drr,     "physical Drr[r^2]=2 test failed");
    }

    // -------------------------------------------------------------------------
    // Test 2: interpolation + physical derivative at off-grid points
    // -------------------------------------------------------------------------

    void test_physical_interpolation()
    {
        std::cout << "[2/3] testing off-grid interpolation and dr-evaluation ...\n";

        constexpr std::size_t Nz = 7;
        constexpr std::size_t Nr = 25;

        SpectralDiffer differ(Nz, Nr);
        ghz::numeric::Domain domain{2.0, 6.0};
        ghz::numeric::PhysicalChebRadialOps rops(differ, domain);

        const auto& r = rops.r();
        const std::size_t N = r.size();

        // Complex quadratic polynomial in physical r:
        // f(r) = c0 + c1 r + c2 r^2
        const Complex c0(1.0,  2.0);
        const Complex c1(3.0, -1.0);
        const Complex c2(0.5,  0.25);

        CVector vals(N);
        for (std::size_t i = 0; i < N; ++i) {
            vals[i] = c0 + c1 * r[i] + c2 * r[i] * r[i];
        }

        const std::vector<Real> test_points{
                Real(2.10), Real(2.75), Real(3.333), Real(4.20), Real(5.90)
        };

        const Real tol = Real(1e-11);

        for (std::size_t k = 0; k < test_points.size(); ++k) {
            const Real r0 = test_points[k];

            const Complex exact_f  = c0 + c1 * r0 + c2 * r0 * r0;
            const Complex exact_df = c1 + Real(2) * c2 * r0;

            const auto out = rops.interp_value_and_dr_from_complex_values(vals, r0);
            const Complex got_f  = out.first;
            const Complex got_df = out.second;

            {
                std::ostringstream os;
                os << "interp f(r) at point " << k;
                require_close_complex(got_f, exact_f, tol, os.str());
            }
            {
                std::ostringstream os;
                os << "interp f'(r) at point " << k;
                require_close_complex(got_df, exact_df, tol, os.str());
            }
        }

        std::cout << "    interpolation and dr evaluation passed.\n";
    }


    inline Real max_residual_inf_norm(
            const matrix<Real>& B,
            const std::vector<Complex>& u,
            const std::vector<Complex>& rhs)
    {
        const std::size_t N = u.size();

        if (B.size1() != N || B.size2() != N || rhs.size() != N) {
            throw std::runtime_error("max_residual_inf_norm: size mismatch.");
        }

        Real max_res = Real(0);

        for (std::size_t i = 0; i < N; ++i) {
            Complex sum = Complex(0.0, 0.0);
            for (std::size_t j = 0; j < N; ++j) {
                sum += B(i, j) * u[j];
            }
            max_res = std::max(max_res, std::abs(sum - rhs[i]));
        }

        return max_res;
    }


    inline Real max_rhs_inf_norm(const std::vector<Complex>& rhs)
    {
        Real max_rhs = Real(0);
        for (const auto& v : rhs) {
            max_rhs = std::max(max_rhs, std::abs(v));
        }
        return max_rhs;
    }

    // -------------------------------------------------------------------------
    // Test 3: manufactured one-domain Schwarzschild collocation solve
    //
    // Solve:
    //    -(u'' + 2/r u') = f
    // with exact solution u(r) = r^3
    //
    // Then:
    //    u'(r) = 3 r^2
    //    u''(r)= 6 r
    //    f(r)  = -(6r + 2/r * 3r^2) = -12 r
    //
    // Left BCs:
    //    u(rmin)  = rmin^3
    //    u'(rmin) = 3 rmin^2
    // -------------------------------------------------------------------------

    void test_schwarzschild_collocation_manufactured()
    {
        std::cout << "[3/3] testing manufactured Schwarzschild collocation solve ...\n";

        constexpr std::size_t Nz = 5;
        constexpr std::size_t Nr = 41;

        SpectralDiffer differ(Nz, Nr);
        ghz::numeric::Domain domain{2.0, 6.0};
        ghz::numeric::PhysicalChebRadialOps rops(differ, domain);

        const auto& r   = rops.r();
        const auto& Dr  = rops.Dr_matrix();
        const auto& Drr = rops.Drr_matrix();

        const std::size_t N = r.size();
        require_true(N >= 3, "Need at least 3 radial points.");

        // Build operator A = -(Drr + 2 diag(1/r) Dr)
        matrix<Real> A(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                A(i, j) = -(Drr(i, j) + Real(2) / r[i] * Dr(i, j));
            }
        }

        // collocation system B u = rhs:
        // row 0: u(rmin) = T1
        // row 1: u'(rmin) = T2
        // rows 2..: A u = f
        matrix<Real> B(A);

        for (std::size_t j = 0; j < N; ++j) {
            B(0, j) = Real(0);
        }
        B(0, 0) = Real(1);

        for (std::size_t j = 0; j < N; ++j) {
            B(1, j) = Dr(0, j);
        }

        std::vector<Complex> rhs(N, Complex(0.0, 0.0));

        const Real rmin = r.front();

        // exact solution u = r^3
        const auto u_exact  = [](Real rr) -> Complex { return Complex(rr * rr * rr, 0.0); };
        const auto up_exact = [](Real rr) -> Complex { return Complex(Real(3) * rr * rr, 0.0); };
        const auto f_exact  = [](Real rr) -> Complex { return Complex(Real(-12) * rr, 0.0); };

        rhs[0] = u_exact(rmin);
        rhs[1] = up_exact(rmin);
        for (std::size_t i = 2; i < N; ++i) {
            rhs[i] = f_exact(r[i]);
        }

        const auto u_num = solve_real_matrix_complex_rhs(B, rhs);
        const Real max_res = max_residual_inf_norm(B, u_num, rhs);

        std::cout << "    max ||B u - rhs||_inf          = "
                  << std::setprecision(18) << max_res << "\n";

        require_true(max_res < Real(1e-10),
                     "manufactured Schwarzschild collocation residual too large");
        Real max_err = Real(0);
        for (std::size_t i = 0; i < N; ++i) {
            max_err = std::max(max_err, std::abs(u_num[i] - u_exact(r[i])));
        }
        const Real rhs_inf = max_rhs_inf_norm(rhs);
        const Real rel_res = (rhs_inf > Real(0)) ? max_res / rhs_inf : max_res;

        std::cout << "    relative residual              = "
                  << std::setprecision(18) << rel_res << "\n";

        std::cout << "    max nodal error in manufactured solution = "
                  << std::setprecision(18) << max_err << "\n";

        require_true(max_err < Real(5e-8),
                     "manufactured Schwarzschild collocation solve failed");
    }

} // namespace

int main()
{
    try {
        test_physical_derivatives();
        test_physical_interpolation();
        test_schwarzschild_collocation_manufactured();

        std::cout << "\nAll spectral physical layer tests passed.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nTest failure: " << e.what() << "\n";
        return 1;
    }
}