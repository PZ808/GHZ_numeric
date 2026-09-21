//
// Created by Peter Zimmerman on 31.10.25.
//
#ifndef GHZ_NUMERIC_CORRECTOR_HPP
#define GHZ_NUMERIC_CORRECTOR_HPP
#pragma once

#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/transport/RK/ODE.hpp"
#include "ghz/transport/collocation/SpectralSliceSolve.hpp"
#include "ghz/geom/DataDomain.hpp"

#include <functional>
#include <utility>
#include <vector>
#include <string>

// structure for solving the transport equations is Sequential+Domainwise
//===================================================
// Sequential in hierarchy
//===================================================
// These quantities depend on earlier solved fields:
// source_xnm depends on Xmmbar, source_xnn depends on Xmmbar and Xnm
//===================================================
// Domain-wise in storage and operators
//===================================================
// Everything is computed separately on .left and .right:
// - edthH / edthBarH / thornPH = Dt operators
// - source assembly
// - collocation solve inputs

namespace ghz::transport {

    using namespace teuk::literals;
    using StateVec =  ode::StateVec;
    using StateMat =  std::vector<std::vector<StateVec>>;
    using RVector = std::vector<Real>;
    using CVector = std::vector<Complex>;
    using GHPSpectral = spectral::SpectralGHPVectorized;
    using Modes = spectral::SpectralGHPVectorized::Modes;
    using GHPType = ghp::GHPType;


    // Layout for the 4-real state [ReX, ReX', ImX, ImX']
    struct Layout4 {
        int re = 0;
        int reR = 1;
        int im  = 2;
        int imR = 3;
    };


    struct TwoDomainField {
        spectral::SpectralGHPVectorized left;
        spectral::SpectralGHPVectorized right;
    };


    template <typename DerivPack>
    struct TwoDomainPack {
        DerivPack left;
        DerivPack right;
    };


    inline TwoDomainField make_two_domain_field(
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

    struct WorldtubeHierarchyResult {
        TwoDomainField xmmbar;
        TwoDomainField xnm;
        TwoDomainField xnn;
    };

    enum class TransportOrder {
        First,
        Second
    };

    // This is the orbit-independent contract for one hierarchy level.  An
    // orbit/source implementation supplies the conditioned equation and
    // boundary data; the solver only sees radial slice data.
    struct WorldtubeLevelCallbacks {
        TransportOrder order{TransportOrder::Second};
        GHPType output_type{};

        using EquationBuilder = std::function<
                collocation::TwoDomainEquationSlice(
                        size_t iz,
                        const TwoDomainField* lower,
                        const TwoDomainField* middle)>;

        EquationBuilder equation;
        // Optional one-time hook after lower hierarchy fields are available.
        // A generic-orbit source adapter can use this to build N[X_mmbar] or
        // U[X_mmbar], V[X_nm] on the complete z/r grid before slicing.
        std::function<void(const TwoDomainField* lower,
                           const TwoDomainField* middle)> prepare;
        std::function<Complex(size_t iz)> bc_value_left;
        std::function<Complex(size_t iz)> bc_deriv_left;
        collocation::InterfaceConditionBuilder interface_conditions;
    };

    // Solve one level over all z slices.  For X_nm, lower is X_mmbar.  For
    // X_nn, lower is X_mmbar and middle is X_nm.
    TwoDomainField solve_worldtube_level(
            const spectral::SpectralDiffer& differ_left,
            const spectral::SpectralDiffer& differ_right,
            const numeric::PunctureTwoDomainSplit& split,
            const RVector& z_grid,
            const Modes& modes,
            const WorldtubeLevelCallbacks& callbacks,
            const TwoDomainField* lower = nullptr,
            const TwoDomainField* middle = nullptr);

    // Sequentially solve X_mmbar, X_nm, and X_nn.  No orbit, archive, or
    // circular-orbit assumptions enter this routine.
    WorldtubeHierarchyResult solve_worldtube_hierarchy(
            const spectral::SpectralDiffer& differ_left,
            const spectral::SpectralDiffer& differ_right,
            const numeric::PunctureTwoDomainSplit& split,
            const RVector& z_grid,
            const Modes& modes,
            const WorldtubeLevelCallbacks& xmmbar,
            const WorldtubeLevelCallbacks& xnm,
            const WorldtubeLevelCallbacks& xnn);


} // namespace ghz::transport


#endif //GHZ_NUMERIC_CORRECTOR_HPP//
