#ifndef GHZ_TRANSPORT_COLLOCATION_SPECTRALSLICESOLVE_HPP
#define GHZ_TRANSPORT_COLLOCATION_SPECTRALSLICESOLVE_HPP

#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/spectral/PhysicalChebRadialOps.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"

#include <boost/numeric/ublas/matrix.hpp>

#include <functional>
#include <vector>

namespace ghz::collocation {

    using teuk::Real;
    using teuk::Complex;

    using SpectralField = spectral::SpectralGHPVectorized;
    using Modes = SpectralField::Modes;

    namespace ublas = boost::numeric::ublas;

    enum class BCKind {
        Value,
        Derivative
    };

    enum class BCSide {
        Left,
        Right
    };

    struct BoundaryCondition {
        BCKind kind{BCKind::Value};
        BCSide side{BCSide::Left};
        Complex value{0.0, 0.0};
    };

    struct ScalarEquationSlice {
        std::vector<Complex> a2;
        std::vector<Complex> a1;
        std::vector<Complex> a0;
        std::vector<Complex> rhs;
    };

    struct TwoDomainEquationSlice {
        ScalarEquationSlice left;
        ScalarEquationSlice right;
    };

    struct InterfaceCondition {
        Complex value_jump{0.0, 0.0};
        Complex derivative_jump{0.0, 0.0};
    };

    struct TwoDomainSolutionSlice {
        std::vector<Complex> left;
        std::vector<Complex> right;
    };

    struct TwoDomainSolutionGrid {
        SpectralField left;
        SpectralField right;
    };

    using TwoDomainSliceEquationBuilder =
            std::function<TwoDomainEquationSlice(size_t iz)>;

    using SliceBCBuilder =
            std::function<BoundaryCondition(size_t iz)>;

    using InterfaceConditionBuilder =
            std::function<InterfaceCondition(size_t iz)>;

    TwoDomainSolutionSlice solve_scalar_two_domain_slice(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            const TwoDomainEquationSlice& eq,
            const BoundaryCondition& bc_row0,
            const BoundaryCondition& bc_row1,
            const InterfaceCondition& iface);

    TwoDomainSolutionGrid solve_scalar_two_domain_grid(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            size_t Nz,
            const Modes& modes,
            ghp::GHPType out_type,
            const TwoDomainSliceEquationBuilder& eq_builder,
            const SliceBCBuilder& bc_row0_builder,
            const SliceBCBuilder& bc_row1_builder,
            const InterfaceConditionBuilder& iface_builder);

    Real max_residual_inf_norm(
            const ublas::matrix<Complex>& B,
            const std::vector<Complex>& u,
            const std::vector<Complex>& rhs);

} // namespace ghz::collocation

#endif // GHZ_TRANSPORT_COLLOCATION_SPECTRALSLICESOLVE_HPP