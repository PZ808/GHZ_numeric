//
// Created by Peter Zimmerman on 19.03.26.
//
//
// spectral_solve_unit_tests.cpp
//
// Unit tests for ghz::collocation::solve_scalar_two_domain_slice
//
// Compile example:
//   c++ -O2 -std=c++20 -fopenmp \
//      spectral_solve_unit_tests.cpp \
//      SpectralDiffer.cpp PhysicalChebRadialOps.cpp SpectralSliceSolve.cpp \
//      -I/path/to/include -o spectral_solve_unit_tests
//

#include "ghz/core/GhzTypes.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/spectral/PhysicalChebRadialOps.hpp"
#include "ghz/transport/collocation/SpectralSliceSolve.hpp"
#include "ghz/geom/DataDomain.hpp"

#include <boost/numeric/ublas/matrix.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <exception>
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
    namespace ublas = boost::numeric::ublas;

    using ghz::collocation::BCKind;
    using ghz::collocation::BCSide;
    using ghz::collocation::BoundaryCondition;
    using ghz::collocation::ScalarEquationSlice;
    using ghz::collocation::TwoDomainEquationSlice;
    using ghz::collocation::InterfaceCondition;
    using ghz::collocation::TwoDomainSolutionSlice;

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

    // -------------------------------------------------------------------------
    // Assemble the same block matrix as SpectralSliceSolve.cpp, for residual checks.
    // -------------------------------------------------------------------------
    ublas::matrix<Complex> assemble_block_matrix(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            const TwoDomainEquationSlice& eq,
            const BoundaryCondition& bc_row0,
            const BoundaryCondition& bc_row1,
            const InterfaceCondition& iface,
            std::vector<Complex>& rhs_out)
    {
        const auto& DrL  = rops_left.Dr_matrix();
        const auto& DrrL = rops_left.Drr_matrix();
        const auto& DrR  = rops_right.Dr_matrix();
        const auto& DrrR = rops_right.Drr_matrix();

        const std::size_t NL = rops_left.r().size();
        const std::size_t NR = rops_right.r().size();
        const std::size_t NT = NL + NR;

        ublas::matrix<Complex> B(NT, NT, Complex(0.0, 0.0));
        rhs_out.assign(NT, Complex(0.0, 0.0));

        // Fill left block
        for (std::size_t i = 0; i < NL; ++i) {
            rhs_out[i] = eq.left.rhs[i];
            for (std::size_t j = 0; j < NL; ++j) {
                const Complex Iij = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
                B(i, j) = eq.left.a2[i] * Complex(DrrL(i, j), 0.0)
                          + eq.left.a1[i] * Complex(DrL(i, j), 0.0)
                          + eq.left.a0[i] * Iij;
            }
        }

        // Fill right block
        for (std::size_t i = 0; i < NR; ++i) {
            rhs_out[NL + i] = eq.right.rhs[i];
            for (std::size_t j = 0; j < NR; ++j) {
                const Complex Iij = (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
                B(NL + i, NL + j) = eq.right.a2[i] * Complex(DrrR(i, j), 0.0)
                                    + eq.right.a1[i] * Complex(DrR(i, j), 0.0)
                                    + eq.right.a0[i] * Iij;
            }
        }

        auto overwrite_outer_bc_row = [&](const BoundaryCondition& bc, std::size_t row) {
            for (std::size_t j = 0; j < NT; ++j) {
                B(row, j) = Complex(0.0, 0.0);
            }

            if (bc.side == BCSide::Left) {
                const std::size_t idxL = 0;
                if (bc.kind == BCKind::Value) {
                    B(row, idxL) = Complex(1.0, 0.0);
                } else {
                    for (std::size_t j = 0; j < NL; ++j) {
                        B(row, j) = Complex(DrL(idxL, j), 0.0);
                    }
                }
            } else {
                const std::size_t idxR = NR - 1;
                if (bc.kind == BCKind::Value) {
                    B(row, NL + idxR) = Complex(1.0, 0.0);
                } else {
                    for (std::size_t j = 0; j < NR; ++j) {
                        B(row, NL + j) = Complex(DrR(idxR, j), 0.0);
                    }
                }
            }

            rhs_out[row] = bc.value;
        };

        overwrite_outer_bc_row(bc_row0, 0);
        overwrite_outer_bc_row(bc_row1, 1);

        // Interface value row at row NL-1
        {
            const std::size_t row = NL - 1;
            const std::size_t iL = NL - 1;
            const std::size_t iR = 0;

            for (std::size_t j = 0; j < NT; ++j) {
                B(row, j) = Complex(0.0, 0.0);
            }
            B(row, iL)      = Complex(-1.0, 0.0);
            B(row, NL + iR) = Complex( 1.0, 0.0);
            rhs_out[row] = iface.value_jump;
        }

        // Interface derivative row at row NL
        {
            const std::size_t row = NL;
            const std::size_t iL = NL - 1;
            const std::size_t iR = 0;

            for (std::size_t j = 0; j < NT; ++j) {
                B(row, j) = Complex(0.0, 0.0);
            }
            for (std::size_t j = 0; j < NL; ++j) {
                B(row, j) = Complex(-DrL(iL, j), 0.0);
            }
            for (std::size_t j = 0; j < NR; ++j) {
                B(row, NL + j) = Complex(DrR(iR, j), 0.0);
            }
            rhs_out[row] = iface.derivative_jump;
        }

        return B;
    }

    Real max_nodal_error(const std::vector<Real>& r,
                         const std::vector<Complex>& u,
                         const std::function<Complex(Real)>& exact)
    {
        Real err = Real(0);
        for (std::size_t i = 0; i < r.size(); ++i) {
            err = std::max(err, std::abs(u[i] - exact(r[i])));
        }
        return err;
    }

    void test_two_domain_cubic_zero_jump()
    {
        std::cout << "[1/2] testing two-domain smooth cubic with zero jumps ...\n";

        constexpr std::size_t Nz_dummy = 7;
        constexpr std::size_t Nr_left  = 25;
        constexpr std::size_t Nr_right = 27;

        SpectralDiffer differL(Nz_dummy, Nr_left);
        SpectralDiffer differR(Nz_dummy, Nr_right);

        ghz::numeric::PunctureTwoDomainSplit split(2.0, 3.5, 5.0);

        ghz::numeric::PhysicalChebRadialOps ropsL(differL, split.left.interval);
        ghz::numeric::PhysicalChebRadialOps ropsR(differR, split.right.interval);

        const auto& rL = ropsL.r();
        const auto& rR = ropsR.r();

        auto exact_u   = [](Real r) -> Complex { return Complex(r*r*r, 0.0); };
        auto exact_du  = [](Real r) -> Complex { return Complex(3.0*r*r, 0.0); };
        auto exact_rhs = [](Real r) -> Complex { return Complex(6.0*r, 0.0); };

        TwoDomainEquationSlice eq;
        eq.left.a2.resize(rL.size(), Complex(1.0, 0.0));
        eq.left.a1.resize(rL.size(), Complex(0.0, 0.0));
        eq.left.a0.resize(rL.size(), Complex(0.0, 0.0));
        eq.left.rhs.resize(rL.size());

        eq.right.a2.resize(rR.size(), Complex(1.0, 0.0));
        eq.right.a1.resize(rR.size(), Complex(0.0, 0.0));
        eq.right.a0.resize(rR.size(), Complex(0.0, 0.0));
        eq.right.rhs.resize(rR.size());

        for (std::size_t i = 0; i < rL.size(); ++i) eq.left.rhs[i] = exact_rhs(rL[i]);
        for (std::size_t i = 0; i < rR.size(); ++i) eq.right.rhs[i] = exact_rhs(rR[i]);

        BoundaryCondition bc0{BCKind::Value,      BCSide::Left, exact_u(rL.front())};
        BoundaryCondition bc1{BCKind::Derivative, BCSide::Left, exact_du(rL.front())};
        InterfaceCondition iface{Complex(0.0, 0.0), Complex(0.0, 0.0)};

        const auto sol = ghz::collocation::solve_scalar_two_domain_slice(
                ropsL, ropsR, eq, bc0, bc1, iface);

        const Real errL = max_nodal_error(rL, sol.left, exact_u);
        const Real errR = max_nodal_error(rR, sol.right, exact_u);

        const auto& DrL = ropsL.Dr_matrix();
        const auto& DrR = ropsR.Dr_matrix();

        Complex u_jump  = sol.right.front() - sol.left.back();

        Complex duL = Complex(0.0, 0.0);
        for (std::size_t j = 0; j < sol.left.size(); ++j) {
            duL += DrL(sol.left.size()-1, j) * sol.left[j];
        }

        Complex duR = Complex(0.0, 0.0);
        for (std::size_t j = 0; j < sol.right.size(); ++j) {
            duR += DrR(0, j) * sol.right[j];
        }

        Complex du_jump = duR - duL;

        std::vector<Complex> rhs_block;
        auto B = assemble_block_matrix(ropsL, ropsR, eq, bc0, bc1, iface, rhs_block);

        std::vector<Complex> u_block(sol.left.size() + sol.right.size());
        for (std::size_t i = 0; i < sol.left.size(); ++i) {
            u_block[i] = sol.left[i];
        }
        for (std::size_t i = 0; i < sol.right.size(); ++i) {
            u_block[sol.left.size() + i] = sol.right[i];
        }

        const Real res_inf = ghz::collocation::max_residual_inf_norm(B, u_block, rhs_block);

        std::cout << "    max left nodal error      = " << std::setprecision(18) << errL << "\n";
        std::cout << "    max right nodal error     = " << std::setprecision(18) << errR << "\n";
        std::cout << "    interface value jump      = " << std::setprecision(18) << std::abs(u_jump) << "\n";
        std::cout << "    interface derivative jump = " << std::setprecision(18) << std::abs(du_jump) << "\n";
        std::cout << "    residual inf norm         = " << std::setprecision(18) << res_inf << "\n";

        require_true(errL < Real(1e-8), "smooth cubic: left nodal error too large");
        require_true(errR < Real(1e-8), "smooth cubic: right nodal error too large");
        require_true(std::abs(u_jump)  < Real(1e-9), "smooth cubic: interface value jump too large");
        require_true(std::abs(du_jump) < Real(1e-8), "smooth cubic: interface derivative jump too large");
        require_true(res_inf < Real(1e-8), "smooth cubic: residual too large");
    }

    void test_two_domain_affine_nonzero_jump()
    {
        std::cout << "[2/2] testing two-domain affine solution with prescribed jumps ...\n";

        constexpr std::size_t Nz_dummy = 7;
        constexpr std::size_t Nr_left  = 21;
        constexpr std::size_t Nr_right = 23;

        SpectralDiffer differL(Nz_dummy, Nr_left);
        SpectralDiffer differR(Nz_dummy, Nr_right);

        ghz::numeric::PunctureTwoDomainSplit split(2.0, 3.5, 5.0);

        ghz::numeric::PhysicalChebRadialOps ropsL(differL, split.left.interval);
        ghz::numeric::PhysicalChebRadialOps ropsR(differR, split.right.interval);

        const auto& rL = ropsL.r();
        const auto& rR = ropsR.r();
        const Real rp = split.r_p;

        const Complex A(2.0, -0.5);
        const Complex B(1.25, 0.75);
        const Complex J0(0.3, -0.2);
        const Complex J1(-0.4, 0.1);

        auto exact_uL  = [&](Real r) -> Complex { return A + B*r; };
        auto exact_duL = [&](Real /*r*/) -> Complex { return B; };

        auto exact_uR  = [&](Real r) -> Complex {
            return exact_uL(r) + J0 + J1*(r - rp);
        };
        auto exact_duR = [&](Real /*r*/) -> Complex {
            return B + J1;
        };

        TwoDomainEquationSlice eq;
        eq.left.a2.resize(rL.size(), Complex(1.0, 0.0));
        eq.left.a1.resize(rL.size(), Complex(0.0, 0.0));
        eq.left.a0.resize(rL.size(), Complex(0.0, 0.0));
        eq.left.rhs.resize(rL.size(), Complex(0.0, 0.0));

        eq.right.a2.resize(rR.size(), Complex(1.0, 0.0));
        eq.right.a1.resize(rR.size(), Complex(0.0, 0.0));
        eq.right.a0.resize(rR.size(), Complex(0.0, 0.0));
        eq.right.rhs.resize(rR.size(), Complex(0.0, 0.0));

        BoundaryCondition bc0{BCKind::Value,      BCSide::Left, exact_uL(rL.front())};
        BoundaryCondition bc1{BCKind::Derivative, BCSide::Left, exact_duL(rL.front())};
        InterfaceCondition iface{J0, J1};

        const auto sol = ghz::collocation::solve_scalar_two_domain_slice(
                ropsL, ropsR, eq, bc0, bc1, iface);

        const Real errL = max_nodal_error(rL, sol.left, exact_uL);
        const Real errR = max_nodal_error(rR, sol.right, exact_uR);

        const auto& DrL = ropsL.Dr_matrix();
        const auto& DrR = ropsR.Dr_matrix();

        Complex u_jump  = sol.right.front() - sol.left.back();

        Complex duL = Complex(0.0, 0.0);
        for (std::size_t j = 0; j < sol.left.size(); ++j) {
            duL += DrL(sol.left.size()-1, j) * sol.left[j];
        }

        Complex duR = Complex(0.0, 0.0);
        for (std::size_t j = 0; j < sol.right.size(); ++j) {
            duR += DrR(0, j) * sol.right[j];
        }

        Complex du_jump = duR - duL;

        std::vector<Complex> rhs_block;
        auto Bmat = assemble_block_matrix(ropsL, ropsR, eq, bc0, bc1, iface, rhs_block);

        std::vector<Complex> u_block(sol.left.size() + sol.right.size());
        for (std::size_t i = 0; i < sol.left.size(); ++i) {
            u_block[i] = sol.left[i];
        }
        for (std::size_t i = 0; i < sol.right.size(); ++i) {
            u_block[sol.left.size() + i] = sol.right[i];
        }

        const Real res_inf = ghz::collocation::max_residual_inf_norm(Bmat, u_block, rhs_block);

        std::cout << "    max left nodal error      = " << std::setprecision(18) << errL << "\n";
        std::cout << "    max right nodal error     = " << std::setprecision(18) << errR << "\n";
        std::cout << "    interface value jump err  = " << std::setprecision(18) << std::abs(u_jump - J0) << "\n";
        std::cout << "    interface deriv jump err  = " << std::setprecision(18) << std::abs(du_jump - J1) << "\n";
        std::cout << "    residual inf norm         = " << std::setprecision(18) << res_inf << "\n";

        require_true(errL < Real(1e-8), "affine jump: left nodal error too large");
        require_true(errR < Real(1e-8), "affine jump: right nodal error too large");
        require_true(std::abs(u_jump  - J0) < Real(1e-9), "affine jump: value jump mismatch");
        require_true(std::abs(du_jump - J1) < Real(1e-8), "affine jump: derivative jump mismatch");
        require_true(res_inf < Real(1e-8), "affine jump: residual too large");
    }

} // namespace

int main()
{
    try {
        test_two_domain_cubic_zero_jump();
        test_two_domain_affine_nonzero_jump();

        std::cout << "\nAll spectral two-domain collocation tests passed.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nTest failure: " << e.what() << "\n";
        return 1;
    }
}