//
// Created by Peter Zimmerman on 19.03.26.
//

//
// xmmbar_worldtube_unit_tests.cpp
//
// Tests the worldtube orchestration helper
//   solve_two_domain_worldtube_grid(...)
// using manufactured two-domain equations.
//

#include "ghz/core/GhzTypes.hpp"
#include "ghz/ghp/GHPScalars.hpp"

#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/spectral/PhysicalChebRadialOps.hpp"
#include "ghz/transport/collocation/SpectralSliceSolve.hpp"
#include "ghz/transport/collocation/XmmbarWorldtubeSolve.hpp"
#include "ghz/geom/DataDomain.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

    using teuk::Real;
    using teuk::Complex;
    using spectral::SpectralDiffer;

    [[noreturn]] void fail(const std::string& msg) {
        throw std::runtime_error(msg);
    }

    void require_true(bool cond, const std::string& msg) {
        if (!cond) fail(msg);
    }

    Real max_nodal_error(
            const std::vector<Real>& r,
            const spectral::SpectralGHPVectorized& field,
            size_t iz,
            const std::function<Complex(Real, Real)>& exact,
            Real z)
    {
        Real err = Real(0);
        for (size_t ir = 0; ir < r.size(); ++ir) {
            err = std::max(err, std::abs(field(ir, iz).value() - exact(r[ir], z)));
        }
        return err;
    }

    void test_worldtube_grid_helper_cubic_family()
    {
        std::cout << "[1/1] testing solve_two_domain_worldtube_grid on manufactured z-dependent cubic family ...\n";

        constexpr size_t Nz = 33;
        constexpr size_t NrL = 65;
        constexpr size_t NrR = 135;

        SpectralDiffer differL(Nz, NrL);
        SpectralDiffer differR(Nz, NrR);

        ghz::numeric::PunctureTwoDomainSplit split(2.0, 3.5, 5.0);

        ghz::numeric::PhysicalChebRadialOps ropsL(differL, split.left.interval);
        ghz::numeric::PhysicalChebRadialOps ropsR(differR, split.right.interval);

        const std::vector<Real> z_grid = differL.lgl_nodes();
        const auto& rL = ropsL.r();
        const auto& rR = ropsR.r();

        ghz::transport::Modes modes{2, 0, 0};
        ghp::GHPType out_type{0, 0};

        // Manufactured family:
        //
        // u(r,z) = alpha(z) * r^3 + beta(z)
        //
        // with
        // alpha(z) = 1 + z
        // beta(z)  = 2 - i z
        //
        // Then
        // u''(r,z) = 6 alpha(z) r
        //
        auto alpha = [](Real z) -> Complex { return Complex(1.0 + z, 0.0); };
        auto beta  = [](Real z) -> Complex { return Complex(2.0, -z); };

        auto exact_u = [&](Real r, Real z) -> Complex {
            return alpha(z) * (r*r*r) + beta(z);
        };

        auto exact_du = [&](Real r, Real z) -> Complex {
            return Complex(3.0 * r * r, 0.0) * alpha(z);
        };

        auto exact_rhs = [&](Real r, Real z) -> Complex {
            return Complex(6.0 * r, 0.0) * alpha(z);
        };

        auto eq_builder = [&](size_t iz) -> ghz::collocation::TwoDomainEquationSlice {
            const Real z = z_grid[iz];

            ghz::collocation::TwoDomainEquationSlice eq;

            eq.left.a2.resize(rL.size(), Complex(1.0, 0.0));
            eq.left.a1.resize(rL.size(), Complex(0.0, 0.0));
            eq.left.a0.resize(rL.size(), Complex(0.0, 0.0));
            eq.left.rhs.resize(rL.size());

            eq.right.a2.resize(rR.size(), Complex(1.0, 0.0));
            eq.right.a1.resize(rR.size(), Complex(0.0, 0.0));
            eq.right.a0.resize(rR.size(), Complex(0.0, 0.0));
            eq.right.rhs.resize(rR.size());

            for (size_t i = 0; i < rL.size(); ++i) {
                eq.left.rhs[i] = exact_rhs(rL[i], z);
            }
            for (size_t i = 0; i < rR.size(); ++i) {
                eq.right.rhs[i] = exact_rhs(rR[i], z);
            }

            return eq;
        };

        auto bc_value_left_builder = [&](size_t iz) -> Complex {
            const Real z = z_grid[iz];
            return exact_u(rL.front(), z);
        };

        auto bc_deriv_left_builder = [&](size_t iz) -> Complex {
            const Real z = z_grid[iz];
            return exact_du(rL.front(), z);
        };

        auto iface_builder = [&](size_t) -> ghz::collocation::InterfaceCondition {
            return ghz::collocation::InterfaceCondition{
                    Complex(0.0, 0.0),
                    Complex(0.0, 0.0)
            };
        };

        const auto out = ghz::transport::solve_two_domain_worldtube_grid(
                ropsL,
                ropsR,
                z_grid,
                modes,
                out_type,
                eq_builder,
                bc_value_left_builder,
                bc_deriv_left_builder,
                iface_builder);

        Real max_err_left = Real(0);
        Real max_err_right = Real(0);

        for (size_t iz = 0; iz < Nz; ++iz) {
            const Real z = z_grid[iz];
            max_err_left = std::max(
                    max_err_left,
                    max_nodal_error(rL, out.left, iz, exact_u, z)
            );
            max_err_right = std::max(
                    max_err_right,
                    max_nodal_error(rR, out.right, iz, exact_u, z)
            );
        }

        std::cout << "    max left nodal error  = " << std::setprecision(18) << max_err_left << "\n";
        std::cout << "    max right nodal error = " << std::setprecision(18) << max_err_right << "\n";

        require_true(max_err_left  < Real(1e-9),
                     "worldtube grid helper: left error too large");
        require_true(max_err_right < Real(1e-9),
                     "worldtube grid helper: right error too large");
    }

} // namespace

int main()
{
    try {
        test_worldtube_grid_helper_cubic_family();
        std::cout << "\nAll worldtube orchestration tests passed.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nTest failure: " << e.what() << "\n";
        return 1;
    }
}