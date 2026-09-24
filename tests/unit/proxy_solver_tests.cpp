#include "ghz/transport/Corrector.hpp"
#include "ghz/transport/collocation/TransportEqnBuilders.hpp"
#include "ghz/spectral/PhysicalChebRadialOps.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/geom/DataDomain.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

    using teuk::Complex;
    using teuk::Real;

    void require_true(bool condition, const char* message)
    {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    template <typename Field, typename Function>
    Real max_error(const Field& field,
                   const std::vector<Real>& r,
                   const std::vector<Real>& z,
                   Function exact)
    {
        Real error = Real(0);
        for (size_t iz = 0; iz < z.size(); ++iz) {
            for (size_t ir = 0; ir < r.size(); ++ir) {
                error = std::max(error,
                                 std::abs(field(ir, iz).value() - exact(r[ir], z[iz])));
            }
        }
        return error;
    }

    void test_proxy_data_full_hierarchy()
    {
        constexpr size_t Nz = 7;
        constexpr size_t NrL = 21;
        constexpr size_t NrR = 25;
        constexpr Real a = Real(0.7);
        constexpr Real k_nm = Real(0.35);
        constexpr Real k_nn_mm = Real(0.2);
        constexpr Real k_nn_nm = Real(-0.15);

        spectral::SpectralDiffer differL(Nz, NrL);
        spectral::SpectralDiffer differR(Nz, NrR);
        ghz::numeric::PunctureTwoDomainSplit split(2.0, 3.5, 5.0);
        ghz::numeric::PhysicalChebRadialOps ropsL(differL, split.left.interval);
        ghz::numeric::PhysicalChebRadialOps ropsR(differR, split.right.interval);
        const std::vector<Real> z_grid = differL.lgl_nodes();
        const ghz::transport::Modes modes{3, 2, 1};
        const ghp::GHPType type{0, 0};

        const auto xmm = [](Real r, Real z) {
            return Complex(r * r * r + z, 0.0);
        };
        const auto xmm_r = [](Real r, Real) {
            return Complex(3.0 * r * r, 0.0);
        };
        const auto xmm_rr = [](Real r, Real) {
            return Complex(6.0 * r, 0.0);
        };
        const auto xnm = [](Real r, Real z) {
            return Complex(r * r + Real(0.3) * z, 0.0);
        };
        const auto xnm_r = [](Real r, Real) {
            return Complex(2.0 * r, 0.0);
        };
        const auto xnm_rr = [](Real, Real) {
            return Complex(2.0, 0.0);
        };
        const auto xnn = [](Real r, Real z) {
            return Complex(r + Real(0.2) * z, 0.0);
        };
        const auto xnn_r = [](Real, Real) {
            return Complex(1.0, 0.0);
        };

        bool xmm_prepared = false;
        bool xnm_prepared = false;
        bool xnn_prepared = false;

        ghz::transport::WorldtubeLevelCallbacks xmm_cb;
        xmm_cb.output_type = type;
        xmm_cb.prepare = [&](const auto* lower, const auto* middle) {
            require_true(lower == nullptr && middle == nullptr,
                         "proxy X_mmbar received dependencies");
            xmm_prepared = true;
        };
        xmm_cb.equation = [&](size_t iz, const auto* lower, const auto* middle) {
            require_true(lower == nullptr && middle == nullptr,
                         "proxy X_mmbar equation received dependencies");
            const std::vector<Complex> zeroL(ropsL.r().size(), Complex(0.0, 0.0));
            const std::vector<Complex> zeroR(ropsR.r().size(), Complex(0.0, 0.0));
            auto eq = ghz::collocation::build_xmmbar_two_domain_eq_slice(
                    ropsL, ropsR, z_grid[iz], a, zeroL, zeroR);
            for (size_t i = 0; i < ropsL.r().size(); ++i) {
                eq.left.rhs[i] = eq.left.a2[i] * xmm_rr(ropsL.r()[i], z_grid[iz])
                                 + eq.left.a1[i] * xmm_r(ropsL.r()[i], z_grid[iz])
                                 + eq.left.a0[i] * xmm(ropsL.r()[i], z_grid[iz]);
            }
            for (size_t i = 0; i < ropsR.r().size(); ++i) {
                eq.right.rhs[i] = eq.right.a2[i] * xmm_rr(ropsR.r()[i], z_grid[iz])
                                  + eq.right.a1[i] * xmm_r(ropsR.r()[i], z_grid[iz])
                                  + eq.right.a0[i] * xmm(ropsR.r()[i], z_grid[iz]);
            }
            return eq;
        };
        xmm_cb.bc_value_left = [&](size_t iz) {
            return xmm(ropsL.r().front(), z_grid[iz]);
        };
        xmm_cb.bc_deriv_left = [&](size_t iz) {
            return xmm_r(ropsL.r().front(), z_grid[iz]);
        };
        xmm_cb.interface_conditions = [](size_t) {
            return ghz::collocation::InterfaceCondition{};
        };

        ghz::transport::WorldtubeLevelCallbacks xnm_cb;
        xnm_cb.output_type = type;
        xnm_cb.prepare = [&](const auto* lower, const auto* middle) {
            require_true(lower != nullptr && middle == nullptr,
                         "proxy X_nm preparation dependencies are wrong");
            xnm_prepared = true;
        };
        xnm_cb.equation = [&](size_t iz, const auto* lower, const auto* middle) {
            require_true(lower != nullptr && middle == nullptr,
                         "proxy X_nm equation dependencies are wrong");
            const std::vector<Complex> zeroL(ropsL.r().size(), Complex(0.0, 0.0));
            const std::vector<Complex> zeroR(ropsR.r().size(), Complex(0.0, 0.0));
            auto eq = ghz::collocation::build_xnm_two_domain_eq_slice(
                    ropsL, ropsR, z_grid[iz], a, zeroL, zeroR);
            for (size_t i = 0; i < ropsL.r().size(); ++i) {
                const Real r = ropsL.r()[i];
                const Real z = z_grid[iz];
                const Complex base = eq.left.a2[i] * xnm_rr(r, z)
                                     + eq.left.a1[i] * xnm_r(r, z)
                                     + eq.left.a0[i] * xnm(r, z);
                eq.left.rhs[i] = base + k_nm *
                                        (lower->left(i, iz).value() - xmm(r, z));
            }
            for (size_t i = 0; i < ropsR.r().size(); ++i) {
                const Real r = ropsR.r()[i];
                const Real z = z_grid[iz];
                const Complex base = eq.right.a2[i] * xnm_rr(r, z)
                                     + eq.right.a1[i] * xnm_r(r, z)
                                     + eq.right.a0[i] * xnm(r, z);
                eq.right.rhs[i] = base + k_nm *
                                         (lower->right(i, iz).value() - xmm(r, z));
            }
            return eq;
        };
        xnm_cb.bc_value_left = [&](size_t iz) {
            return xnm(ropsL.r().front(), z_grid[iz]);
        };
        xnm_cb.bc_deriv_left = [&](size_t iz) {
            return xnm_r(ropsL.r().front(), z_grid[iz]);
        };
        xnm_cb.interface_conditions = [](size_t) {
            return ghz::collocation::InterfaceCondition{};
        };

        ghz::transport::WorldtubeLevelCallbacks xnn_cb;
        xnn_cb.order = ghz::transport::TransportOrder::First;
        xnn_cb.output_type = type;
        xnn_cb.prepare = [&](const auto* lower, const auto* middle) {
            require_true(lower != nullptr && middle != nullptr,
                         "proxy X_nn preparation dependencies are wrong");
            xnn_prepared = true;
        };
        xnn_cb.equation = [&](size_t iz, const auto* lower, const auto* middle) {
            require_true(lower != nullptr && middle != nullptr,
                         "proxy X_nn equation dependencies are wrong");
            const std::vector<Complex> zeroL(ropsL.r().size(), Complex(0.0, 0.0));
            const std::vector<Complex> zeroR(ropsR.r().size(), Complex(0.0, 0.0));
            auto eq = ghz::collocation::build_xnn_two_domain_eq_slice(
                    ropsL, ropsR, z_grid[iz], a, zeroL, zeroR);
            for (size_t i = 0; i < ropsL.r().size(); ++i) {
                const Real r = ropsL.r()[i];
                const Real z = z_grid[iz];
                const Complex base = eq.left.a1[i] * xnn_r(r, z)
                                     + eq.left.a0[i] * xnn(r, z);
                eq.left.rhs[i] = base
                                 + k_nn_mm * (lower->left(i, iz).value() - xmm(r, z))
                                 + k_nn_nm * (middle->left(i, iz).value() - xnm(r, z));
            }
            for (size_t i = 0; i < ropsR.r().size(); ++i) {
                const Real r = ropsR.r()[i];
                const Real z = z_grid[iz];
                const Complex base = eq.right.a1[i] * xnn_r(r, z)
                                     + eq.right.a0[i] * xnn(r, z);
                eq.right.rhs[i] = base
                                  + k_nn_mm * (lower->right(i, iz).value() - xmm(r, z))
                                  + k_nn_nm * (middle->right(i, iz).value() - xnm(r, z));
            }
            return eq;
        };
        xnn_cb.bc_value_left = [&](size_t iz) {
            return xnn(ropsL.r().front(), z_grid[iz]);
        };
        xnn_cb.interface_conditions = [](size_t) {
            return ghz::collocation::InterfaceCondition{};
        };

        const auto result = ghz::transport::solve_worldtube_hierarchy(
                differL, differR, split, z_grid, modes,
                xmm_cb, xnm_cb, xnn_cb);

        const Real error = std::max({
                max_error(result.xmmbar.left, ropsL.r(), z_grid, xmm),
                max_error(result.xmmbar.right, ropsR.r(), z_grid, xmm),
                max_error(result.xnm.left, ropsL.r(), z_grid, xnm),
                max_error(result.xnm.right, ropsR.r(), z_grid, xnm),
                max_error(result.xnn.left, ropsL.r(), z_grid, xnn),
                max_error(result.xnn.right, ropsR.r(), z_grid, xnn)});

        require_true(xmm_prepared && xnm_prepared && xnn_prepared,
                     "proxy hierarchy preparation hooks were not called");
        require_true(error < Real(1e-7),
                     "proxy-data hierarchy error is too large");
        std::cout << "proxy hierarchy maximum error = "
                  << std::setprecision(18) << error << "\n";
    }

} // namespace

int main()
{
    try {
        test_proxy_data_full_hierarchy();
        std::cout << "Proxy solver test passed.\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "Proxy solver test FAILED: " << error.what() << '\n';
        return 1;
    }
}
