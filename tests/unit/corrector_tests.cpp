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
#include <functional>
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
        std::cout << "[1/3] testing solve_two_domain_worldtube_grid on manufactured z-dependent cubic family ...\n";

        constexpr size_t Nz = 9;
        constexpr size_t NrL = 25;
        constexpr size_t NrR = 29;

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

    void test_complete_hierarchy_manufactured()
    {
        std::cout << "[2/3] testing sequential X_mmbar -> X_nm -> X_nn hierarchy ...\n";

        constexpr size_t Nz = 17;
        constexpr size_t NrL = 33;
        constexpr size_t NrR = 41;

        SpectralDiffer differL(Nz, NrL);
        SpectralDiffer differR(Nz, NrR);
        ghz::numeric::PunctureTwoDomainSplit split(2.0, 3.5, 5.0);
        ghz::transport::Modes modes{2, 0, 0};
        ghp::GHPType type{0, 0};
        const auto z_grid = differL.lgl_nodes();
        ghz::numeric::PhysicalChebRadialOps ropsL(differL, split.left.interval);
        ghz::numeric::PhysicalChebRadialOps ropsR(differR, split.right.interval);

        const auto exact_xmmbar = [](Real r, Real z) -> Complex {
            return Complex(r * r * r + z, 0.0);
        };
        const auto exact_xnm = [](Real r, Real z) -> Complex {
            return Complex(r * r + Real(0.25) * z, 0.0);
        };
        const auto exact_xnn = [](Real r, Real z) -> Complex {
            return Complex(r + z, 0.0);
        };

        const auto make_second_order =
                [&](auto exact, auto rhs, auto deriv) {
                    ghz::transport::WorldtubeLevelCallbacks cb;
                    cb.order = ghz::transport::TransportOrder::Second;
                    cb.output_type = type;
                    cb.equation = [&, exact, rhs](size_t iz,
                                      const ghz::transport::TwoDomainField* lower,
                                      const ghz::transport::TwoDomainField* middle) {
                        require_true(lower == nullptr && middle == nullptr,
                                     "seed level received unexpected dependencies");
                        ghz::collocation::TwoDomainEquationSlice eq;
                        eq.left.a2.assign(ropsL.r().size(), Complex(1.0, 0.0));
                        eq.left.a1.assign(ropsL.r().size(), Complex(0.0, 0.0));
                        eq.left.a0.assign(ropsL.r().size(), Complex(0.0, 0.0));
                        eq.right.a2.assign(ropsR.r().size(), Complex(1.0, 0.0));
                        eq.right.a1.assign(ropsR.r().size(), Complex(0.0, 0.0));
                        eq.right.a0.assign(ropsR.r().size(), Complex(0.0, 0.0));
                        eq.left.rhs.resize(ropsL.r().size());
                        eq.right.rhs.resize(ropsR.r().size());
                        for (size_t i = 0; i < ropsL.r().size(); ++i) {
                            eq.left.rhs[i] = rhs(ropsL.r()[i], z_grid[iz]);
                        }
                        for (size_t i = 0; i < ropsR.r().size(); ++i) {
                            eq.right.rhs[i] = rhs(ropsR.r()[i], z_grid[iz]);
                        }
                        return eq;
                    };
                    cb.bc_value_left = [&, exact](size_t iz) {
                        return exact(ropsL.r().front(), z_grid[iz]);
                    };
                    cb.bc_deriv_left = [&, deriv](size_t iz) {
                        return deriv(ropsL.r().front(), z_grid[iz]);
                    };
                    cb.interface_conditions = [](size_t) {
                        return ghz::collocation::InterfaceCondition{};
                    };
                    return cb;
                };

        auto xmmbar = make_second_order(
                exact_xmmbar,
                [](Real r, Real) { return Complex(6.0 * r, 0.0); },
                [](Real r, Real) { return Complex(3.0 * r * r, 0.0); });

        auto xnm = make_second_order(
                exact_xnm,
                [](Real, Real) { return Complex(2.0, 0.0); },
                [](Real r, Real) { return Complex(2.0 * r, 0.0); });
        xnm.equation = [&](size_t iz,
                           const ghz::transport::TwoDomainField* lower,
                           const ghz::transport::TwoDomainField* middle) {
            require_true(lower != nullptr && middle == nullptr,
                         "X_nm dependency contract was not honored");
            ghz::collocation::TwoDomainEquationSlice eq;
            eq.left.a2.assign(ropsL.r().size(), Complex(1.0, 0.0));
            eq.left.a1.assign(ropsL.r().size(), Complex(0.0, 0.0));
            eq.left.a0.assign(ropsL.r().size(), Complex(0.0, 0.0));
            eq.right.a2.assign(ropsR.r().size(), Complex(1.0, 0.0));
            eq.right.a1.assign(ropsR.r().size(), Complex(0.0, 0.0));
            eq.right.a0.assign(ropsR.r().size(), Complex(0.0, 0.0));
            eq.left.rhs.assign(ropsL.r().size(), Complex(2.0, 0.0));
            eq.right.rhs.assign(ropsR.r().size(), Complex(2.0, 0.0));
            return eq;
        };
        xnm.bc_deriv_left = [&](size_t iz) {
            return Complex(2.0 * ropsL.r().front(), 0.0);
        };

        auto xnn = make_second_order(
                exact_xnn,
                [](Real, Real) { return Complex(1.0, 0.0); },
                [](Real, Real) { return Complex(1.0, 0.0); });
        xnn.order = ghz::transport::TransportOrder::First;
        xnn.bc_deriv_left = {};
        xnn.equation = [&](size_t iz,
                           const ghz::transport::TwoDomainField* lower,
                           const ghz::transport::TwoDomainField* middle) {
            require_true(lower != nullptr && middle != nullptr,
                         "X_nn dependency contract was not honored");
            ghz::collocation::TwoDomainEquationSlice eq;
            eq.left.a2.assign(ropsL.r().size(), Complex(0.0, 0.0));
            eq.left.a1.assign(ropsL.r().size(), Complex(1.0, 0.0));
            eq.left.a0.assign(ropsL.r().size(), Complex(0.0, 0.0));
            eq.right.a2.assign(ropsR.r().size(), Complex(0.0, 0.0));
            eq.right.a1.assign(ropsR.r().size(), Complex(1.0, 0.0));
            eq.right.a0.assign(ropsR.r().size(), Complex(0.0, 0.0));
            eq.left.rhs.assign(ropsL.r().size(), Complex(1.0, 0.0));
            eq.right.rhs.assign(ropsR.r().size(), Complex(1.0, 0.0));
            return eq;
        };

        const auto result = ghz::transport::solve_worldtube_hierarchy(
                differL, differR, split, z_grid, modes, xmmbar, xnm, xnn);

        Real err = Real(0);
        for (size_t iz = 0; iz < z_grid.size(); ++iz) {
            for (size_t ir = 0; ir < ropsL.r().size(); ++ir) {
                err = std::max(err, std::abs(result.xmmbar.left(ir, iz).value()
                                             - exact_xmmbar(ropsL.r()[ir], z_grid[iz])));
                err = std::max(err, std::abs(result.xnm.left(ir, iz).value()
                                             - exact_xnm(ropsL.r()[ir], z_grid[iz])));
                err = std::max(err, std::abs(result.xnn.left(ir, iz).value()
                                             - exact_xnn(ropsL.r()[ir], z_grid[iz])));
            }
        }
        require_true(err < Real(1e-7), "complete hierarchy manufactured error too large");
        std::cout << "    maximum hierarchy error = " << std::setprecision(18) << err << "\n";
    }

    void test_physical_hierarchy_coefficients()
    {
        std::cout << "[3/3] testing physical X_nm and X_nn radial equations ...\n";

        SpectralDiffer differL(9, 25);
        SpectralDiffer differR(9, 29);
        ghz::numeric::PunctureTwoDomainSplit split(2.0, 3.5, 5.0);
        ghz::numeric::PhysicalChebRadialOps ropsL(differL, split.left.interval);
        ghz::numeric::PhysicalChebRadialOps ropsR(differR, split.right.interval);
        const Real z = Real(0.23);
        const Real a = Real(0.7);

        auto rho = [&](Real r) {
            return -Complex(1.0, 0.0) / Complex(r, -a * z);
        };
        auto nm_exact = [&](Real r) {
            const Complex rr = rho(r);
            return rr * (rr + std::conj(rr));
        };
        auto nm_deriv = [&](Real r) {
            const Complex rr = rho(r);
            const Complex rb = std::conj(rr);
            return rr * rr * (rr + rb) + rr * (rr * rr + rb * rb);
        };
        const auto nm_eq = ghz::collocation::build_xnm_two_domain_eq_slice(
                ropsL, ropsR, z, a,
                std::vector<Complex>(ropsL.r().size(), Complex(0.0, 0.0)),
                std::vector<Complex>(ropsR.r().size(), Complex(0.0, 0.0)));
        const auto nm_sol = ghz::collocation::solve_scalar_two_domain_slice(
                ropsL, ropsR, nm_eq,
                {ghz::collocation::BCKind::Value, ghz::collocation::BCSide::Left,
                 nm_exact(ropsL.r().front())},
                {ghz::collocation::BCKind::Derivative, ghz::collocation::BCSide::Left,
                 nm_deriv(ropsL.r().front())},
                {});

        auto nn_exact = [&](Real r) {
            const Complex rr = rho(r);
            return rr + std::conj(rr);
        };
        const auto nn_eq = ghz::collocation::build_xnn_two_domain_eq_slice(
                ropsL, ropsR, z, a,
                std::vector<Complex>(ropsL.r().size(), Complex(0.0, 0.0)),
                std::vector<Complex>(ropsR.r().size(), Complex(0.0, 0.0)));
        const auto nn_sol = ghz::collocation::solve_first_order_two_domain_slice(
                ropsL, ropsR, nn_eq,
                {ghz::collocation::BCKind::Value, ghz::collocation::BCSide::Left,
                 nn_exact(ropsL.r().front())}, {});

        Real max_err = Real(0);
        Real max_nm_err = Real(0);
        Real max_nn_err = Real(0);
        for (size_t i = 0; i < ropsL.r().size(); ++i) {
            max_nm_err = std::max(max_nm_err, std::abs(nm_sol.left[i] - nm_exact(ropsL.r()[i])));
            max_nn_err = std::max(max_nn_err, std::abs(nn_sol.left[i] - nn_exact(ropsL.r()[i])));
        }
        for (size_t i = 0; i < ropsR.r().size(); ++i) {
            max_nm_err = std::max(max_nm_err, std::abs(nm_sol.right[i] - nm_exact(ropsR.r()[i])));
            max_nn_err = std::max(max_nn_err, std::abs(nn_sol.right[i] - nn_exact(ropsR.r()[i])));
        }
        max_err = std::max(max_nm_err, max_nn_err);
        std::cout << "    X_nm error = " << max_nm_err
                  << ", X_nn error = " << max_nn_err << "\n";
        require_true(max_err < Real(1e-7),
                     "physical hierarchy coefficient test error too large");
        std::cout << "    maximum physical-equation error = "
                  << std::setprecision(18) << max_err << "\n";
    }

} // namespace

int main()
{
    try {
        test_worldtube_grid_helper_cubic_family();
        test_complete_hierarchy_manufactured();
        test_physical_hierarchy_coefficients();
        std::cout << "\nAll worldtube orchestration tests passed.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\nTest failure: " << e.what() << "\n";
        return 1;
    }
}
