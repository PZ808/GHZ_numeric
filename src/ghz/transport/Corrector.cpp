// Created by Peter Zimmerman on 01.12.25.
//
#include "ghz/core/GhzTypes.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/transport/Corrector.hpp"

#include <stdexcept>
#include <optional>
#include <string>
namespace ghz::transport {



    using namespace teuk::literals;
    Complex I = teuk::I;

    namespace {

        void validate_level_callbacks(const WorldtubeLevelCallbacks& callbacks,
                                      const char* name)
        {
            if (!callbacks.equation) {
                throw std::invalid_argument(std::string(name) +
                                            ": equation callback is not set");
            }
            if (!callbacks.bc_value_left) {
                throw std::invalid_argument(std::string(name) +
                                            ": left boundary callback is not set");
            }
            if (callbacks.order == TransportOrder::Second &&
                !callbacks.bc_deriv_left) {
                throw std::invalid_argument(std::string(name) +
                                            ": derivative boundary callback is not set");
            }
            if (!callbacks.interface_conditions) {
                throw std::invalid_argument(std::string(name) +
                                            ": interface callback is not set");
            }
        }

        void fill_result_field(
                const std::vector<Complex>& left,
                const std::vector<Complex>& right,
                size_t iz,
                GHPType type,
                TwoDomainField& out)
        {
            for (size_t ir = 0; ir < left.size(); ++ir) {
                out.left.set_index(ir, iz,
                                   ghp::GHPScalar<Complex>(left[ir], type.p, type.q));
            }
            for (size_t ir = 0; ir < right.size(); ++ir) {
                out.right.set_index(ir, iz,
                                    ghp::GHPScalar<Complex>(right[ir], type.p, type.q));
            }
        }

    } // namespace

    TwoDomainField solve_worldtube_level(
            const spectral::SpectralDiffer& differ_left,
            const spectral::SpectralDiffer& differ_right,
            const numeric::PunctureTwoDomainSplit& split,
            const RVector& z_grid,
            const Modes& modes,
            const WorldtubeLevelCallbacks& callbacks,
            const TwoDomainField* lower,
            const TwoDomainField* middle)
    {
        validate_level_callbacks(callbacks, "solve_worldtube_level");
        if (z_grid.empty()) {
            throw std::invalid_argument("solve_worldtube_level: z_grid is empty");
        }

        const numeric::PhysicalChebRadialOps rops_left(
                differ_left, split.left.interval);
        const numeric::PhysicalChebRadialOps rops_right(
                differ_right, split.right.interval);
        TwoDomainField out = make_two_domain_field(
                rops_left.r().size(), rops_right.r().size(), z_grid.size(),
                modes, callbacks.output_type);

        if (callbacks.prepare) {
            callbacks.prepare(lower, middle);
        }

        for (size_t iz = 0; iz < z_grid.size(); ++iz) {
            const auto eq = callbacks.equation(iz, lower, middle);
            const auto iface = callbacks.interface_conditions(iz);
            collocation::TwoDomainSolutionSlice sol;

            if (callbacks.order == TransportOrder::First) {
                sol = collocation::solve_first_order_two_domain_slice(
                        rops_left, rops_right, eq,
                        collocation::BoundaryCondition{
                                collocation::BCKind::Value,
                                collocation::BCSide::Left,
                                callbacks.bc_value_left(iz)},
                        iface);
            } else {
                sol = collocation::solve_scalar_two_domain_slice(
                        rops_left, rops_right, eq,
                        collocation::BoundaryCondition{
                                collocation::BCKind::Value,
                                collocation::BCSide::Left,
                                callbacks.bc_value_left(iz)},
                        collocation::BoundaryCondition{
                                collocation::BCKind::Derivative,
                                collocation::BCSide::Left,
                                callbacks.bc_deriv_left(iz)},
                        iface);
            }

            fill_result_field(sol.left, sol.right, iz,
                              callbacks.output_type, out);
        }
        return out;
    }

    WorldtubeHierarchyResult solve_worldtube_hierarchy(
            const spectral::SpectralDiffer& differ_left,
            const spectral::SpectralDiffer& differ_right,
            const numeric::PunctureTwoDomainSplit& split,
            const RVector& z_grid,
            const Modes& modes,
            const WorldtubeLevelCallbacks& xmmbar,
            const WorldtubeLevelCallbacks& xnm,
            const WorldtubeLevelCallbacks& xnn)
    {
        const numeric::PhysicalChebRadialOps rops_left(
                differ_left, split.left.interval);
        const numeric::PhysicalChebRadialOps rops_right(
                differ_right, split.right.interval);

        WorldtubeHierarchyResult result{
                solve_worldtube_level(differ_left, differ_right, split, z_grid,
                                      modes, xmmbar),
                make_two_domain_field(rops_left.r().size(), rops_right.r().size(),
                                      z_grid.size(), modes, xnm.output_type),
                make_two_domain_field(rops_left.r().size(), rops_right.r().size(),
                                      z_grid.size(), modes, xnn.output_type)};
        result.xnm = solve_worldtube_level(
                differ_left, differ_right, split, z_grid, modes, xnm,
                &result.xmmbar, nullptr);
        result.xnn = solve_worldtube_level(
                differ_left, differ_right, split, z_grid, modes, xnn,
                &result.xmmbar, &result.xnm);
        return result;
    }

} // namespace ghz
