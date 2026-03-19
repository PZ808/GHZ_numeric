//
// Created by Peter Zimmerman on 02.03.26.
//

#ifndef GHZ_NUMERIC_SOURCETORUSFREQUENCIES_HPP
#define GHZ_NUMERIC_SOURCETORUSFREQUENCIES_HPP


// SourceFrequencies.hpp
#pragma once
#include "ghz/core/GhzTypes.hpp"   // Real, Complex, Modes, etc.


namespace orbit {
    using Real = teuk::Real;

    struct SourceFrequencies {
        Real Omega_phi{0};
        Real Omega_r{0};
        Real Omega_z{0};
    };

} //orbit
#endif //GHZ_NUMERIC_SOURCETORUSFREQUENCIES_HPP
