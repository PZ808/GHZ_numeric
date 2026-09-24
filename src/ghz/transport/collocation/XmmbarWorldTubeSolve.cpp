//
// Created by Peter Zimmerman on 19.03.26.
//
#include "ghz/transport/collocation/XmmbarWorldtubeSolve.hpp"

#include <stdexcept>

namespace ghz::transport {

    namespace {

        inline TwoDomainField make_two_domain_field_local(
                size_t Nr_left,
                size_t Nr_right,
                size_t Nz,
                const Modes& modes,
                GHPType type)
        {
            return TwoDomainField{
                    GHPSpectral(
                            Nr_left, Nz, modes,
                            ghp::GHPScalar<Complex>(teuk::zeroC, type.p, type.q),
                            type.p, type.q
                    ),
                    GHPSpectral(
                            Nr_right, Nz, modes,
                            ghp::GHPScalar<Complex>(teuk::zeroC, type.p, type.q),
                            type.p, type.q
                    )
            };
        }

    } // namespace


    TwoDomainField solve_two_domain_worldtube_grid(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            const RVector& z_grid,
            const Modes& modes,
            GHPType out_type,
            const ghz::collocation::TwoDomainSliceEquationBuilder& eq_builder,
            const BoundaryValueBuilder& bc_value_left_builder,
            const BoundaryValueBuilder& bc_deriv_left_builder,
            const ghz::collocation::InterfaceConditionBuilder& iface_builder)
    {

        if (z_grid.empty()) {
            throw std::invalid_argument("solve_two_domain_worldtube_grid: z_grid is empty.");
        }
        if (!eq_builder) {
            throw std::invalid_argument("solve_two_domain_worldtube_grid: eq_builder not set.");
        }
        if (!bc_value_left_builder) {
            throw std::invalid_argument(
                    "solve_two_domain_worldtube_grid: bc_value_left_builder not set.");
        }
        if (!bc_deriv_left_builder) {
            throw std::invalid_argument(
                    "solve_two_domain_worldtube_grid: bc_deriv_left_builder not set.");
        }
        if (!iface_builder) {
            throw std::invalid_argument(
                    "solve_two_domain_worldtube_grid: iface_builder not set.");
        }

        const size_t Nz = z_grid.size();
        const size_t NL = rops_left.r().size();
        const size_t NR = rops_right.r().size();

        TwoDomainField out = make_two_domain_field_local(
                NL, NR, Nz, modes, out_type);

        for (size_t iz = 0; iz < Nz; ++iz) {
            const auto eq = eq_builder(iz);

            const ghz::collocation::BoundaryCondition bc0{
                    ghz::collocation::BCKind::Value,
                    ghz::collocation::BCSide::Left,
                    bc_value_left_builder(iz)
            };

            const ghz::collocation::BoundaryCondition bc1{
                    ghz::collocation::BCKind::Derivative,
                    ghz::collocation::BCSide::Left,
                    bc_deriv_left_builder(iz)
            };

            const auto iface = iface_builder(iz);

            const auto sol = ghz::collocation::solve_scalar_two_domain_slice(
                    rops_left, rops_right, eq, bc0, bc1, iface);

            for (size_t ir = 0; ir < NL; ++ir) {
                out.left.set_index(
                        ir, iz,
                        ghp::GHPScalar<Complex>(sol.left[ir], out_type.p, out_type.q)
                );
            }
            for (size_t ir = 0; ir < NR; ++ir) {
                out.right.set_index(
                        ir, iz,
                        ghp::GHPScalar<Complex>(sol.right[ir], out_type.p, out_type.q)
                );
            }
        }

        return out;
    }

    TwoDomainField solve_xmmbar_worldtube(
            const spectral::SpectralDiffer& differ_left,
            const spectral::SpectralDiffer& differ_right,
            const ghz::numeric::PunctureTwoDomainSplit& split,
            const KerrMetricOutgoing& metric,
            const RVector& z_grid,
            const Modes& modes,
            GHPType out_type,
            ghz::source::ConditionSource& src_left,
            ghz::source::ConditionSource& src_right,
            ghz::source::ConditionBoundaryData& bdy_left,
            const ghz::collocation::InterfaceConditionBuilder& iface_builder)
    {
        if (z_grid.empty()) {
            throw std::invalid_argument("solve_xmmbar_worldtube: z_grid is empty.");
        }
        if (!iface_builder) {
            throw std::invalid_argument("solve_xmmbar_worldtube: iface_builder not set.");
        }

        if (src_left.m() != modes.m) {
            throw std::runtime_error(
                    "solve_xmmbar_worldtube: left source m does not match field modes.m.");
        }
        if (src_right.m() != modes.m) {
            throw std::runtime_error(
                    "solve_xmmbar_worldtube: right source m does not match field modes.m.");
        }
        if (bdy_left.m() != modes.m) {
            throw std::runtime_error(
                    "solve_xmmbar_worldtube: boundary data m does not match field modes.m.");
        }

        ghz::numeric::PhysicalChebRadialOps rops_left(
                differ_left, split.left.interval);
        ghz::numeric::PhysicalChebRadialOps rops_right(
                differ_right, split.right.interval);

        auto eq_builder = [&](size_t iz) -> ghz::collocation::TwoDomainEquationSlice {
            const Real z = z_grid[iz];
            src_left.set_z_slice(z);
            src_right.set_z_slice(z);

            return ghz::collocation::build_xmmbar_two_domain_eq_slice_prepared(
                    rops_left, rops_right, z, metric.a(), src_left, src_right);
        };

        auto bc_value_left_builder = [&](size_t iz) -> Complex {
            const Real z = z_grid[iz];
            bdy_left.set_z_slice(z);
            return bdy_left.delta();
        };

        auto bc_deriv_left_builder = [&](size_t iz) -> Complex {
            const Real z = z_grid[iz];
            bdy_left.set_z_slice(z);
            return bdy_left.deltaPrime();
        };

        return solve_two_domain_worldtube_grid(
                rops_left,
                rops_right,
                z_grid,
                modes,
                out_type,
                eq_builder,
                bc_value_left_builder,
                bc_deriv_left_builder,
                iface_builder);
    }
    //
    // zero jump interface condition overload
    //
    TwoDomainField solve_xmmbar_worldtube(
            const spectral::SpectralDiffer& differ_left,
            const spectral::SpectralDiffer& differ_right,
            const ghz::numeric::PunctureTwoDomainSplit& split,
            const KerrMetricOutgoing& metric,
            const RVector& z_grid,
            const Modes& modes,
            GHPType out_type,
            ghz::source::ConditionSource& src_left,
            ghz::source::ConditionSource& src_right,
            ghz::source::ConditionBoundaryData& bdy_left)
    {
        return solve_xmmbar_worldtube(
                differ_left,
                differ_right,
                split,
                metric,
                z_grid,
                modes,
                out_type,
                src_left,
                src_right,
                bdy_left,
                [](size_t) {
                    return ghz::collocation::InterfaceCondition{
                            Complex(0.0, 0.0),
                            Complex(0.0, 0.0)
                    };
                });
    }

} // namespace ghz::transport
