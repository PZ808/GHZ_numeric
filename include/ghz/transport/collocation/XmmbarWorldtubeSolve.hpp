//
// Created by Peter Zimmerman on 19.03.26.
//
#ifndef GHZ_NUMERIC_TRANSPORT_XMMBARWORLDTUBESOLVE_HPP
#define GHZ_NUMERIC_TRANSPORT_XMMBARWORLDTUBESOLVE_HPP
#pragma once

#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/KerrMetricOutgoing.hpp"
#include "ghz/geom/KerrParams.hpp"
#include "ghz/transport/Corrector.hpp"
#include "ghz/transport/collocation/SpectralSliceSolve.hpp"
#include "ghz/transport/collocation/TransportEqnBuilders.hpp"
#include "ghz/source/ConditionSource.hpp"
#include "ghz/source/ConditionBoundaryData.hpp"
#include "ghz/geom/DataDomain.hpp"

#include <functional>

namespace ghz::transport {

    using BoundaryValueBuilder = std::function<Complex(size_t)>;

    /**
     * @brief Solve X_{m\bar m} on the two puncture domains of the worldtube.
     *
     * Domains:
     *   left  = [r_min, r_p]
     *   right = [r_p,   r_max]
     *
     * For each fixed z_grid[iz], this routine:
     *   1. conditions the left/right effective source on that z-slice,
     *   2. builds the two-domain collocation equation,
     *   3. imposes the two outer boundary rows at r_min,
     *   4. imposes interface jump data at r_p,
     *   5. solves the coupled two-domain block system,
     *   6. stores the result in a TwoDomainField.
     *
     * @note This routine is currently serial in iz because ConditionSource
     *       owns mutable per-slice state via set_z_slice(z).
     */
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
            const ghz::collocation::InterfaceConditionBuilder& iface_builder);

    /**
     * @brief Convenience overload with zero interface jumps.
     */
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
            ghz::source::ConditionBoundaryData& bdy_left);

    TwoDomainField solve_two_domain_worldtube_grid(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            const RVector& z_grid,
            const Modes& modes,
            GHPType out_type,
            const ghz::collocation::TwoDomainSliceEquationBuilder& eq_builder,
            const BoundaryValueBuilder& bc_value_left_builder,
            const BoundaryValueBuilder& bc_deriv_left_builder,
            const ghz::collocation::InterfaceConditionBuilder& iface_builder);

} // namespace ghz::transport

#endif //GHZ_NUMERIC_TRANSPORT_XMMBARWORLDTUBESOLVE_HPP

