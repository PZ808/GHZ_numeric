//
// Created by Peter Zimmerman on 31.10.25.
//


#ifndef GHZ_NUMERIC_CORRECTOR_HPP
#define GHZ_NUMERIC_CORRECTOR_HPP

#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/ghp/GHPFieldVectorized.hpp"
#include "ghz/geom/KerrMetricOutgoing.hpp"
#include "ghz/geom/KinnersleyTetrad.hpp"
#include "ghz/spectral/KinnersleySpectralHeldOperators.hpp"
#include "ODE.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"  // your node builders

#include <utility>
#include <vector>
#include <string>

namespace ghz {

    using namespace teuk::literals;
    using StateVec =  ode::StateVec;
    using StateMat =  std::vector<std::vector<StateVec>>;
    using RVector = std::vector<Real>;
    using CVector = std::vector<Complex>;
    using GHPSpectral = spectral::SpectralGHPVectorized;

    // Layout for the 4-real state [ReX, ReX', ImX, ImX']
    struct Layout4 {
        int re = 0;
        int reR = 1;
        int im  = 2;
        int imR = 3;
    };


} // namespace ghz


#endif //GHZ_NUMERIC_CORRECTOR_HPP//

