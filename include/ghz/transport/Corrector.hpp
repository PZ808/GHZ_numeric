//
// Created by Peter Zimmerman on 31.10.25.
//
#ifndef GHZ_NUMERIC_CORRECTOR_HPP
#define GHZ_NUMERIC_CORRECTOR_HPP
#pragma once

#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/transport/RK/ODE.hpp"

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


} // namespace ghz::transport


#endif //GHZ_NUMERIC_CORRECTOR_HPP//

