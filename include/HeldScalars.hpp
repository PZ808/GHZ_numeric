//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_HELDSCALARS_HPP
#define GHZ_NUMERIC_HELDSCALARS_HPP

#pragma once
#include <complex>
#include <string>
#include "GHPScalars.hpp"

using Complex = std::complex<double>;


class HeldScalarsTypeD {
    // GHP scalars for type-D spacetime in a Weyl aligned NP frame
    // (shear-free geodesic null congruence)
public:
    // Held scalars (unprimed and primed)
    Complex rhopH, tauH,  rhopH_bar, tauH_bar;
    Complex PsiH, PsiH_bar;
    Complex OmH, OmH_bar;


    GHPType rhopH_GHPType{-2,-2};
    GHPType rhopH_bar_GHPType{-2,-2};

    GHPType tauH_GHPType{-1, -3};
    GHPType tauH_bar_GHPType{-3, -1};

    GHPType PsiH_GHPType{-3, -3};
    GHPType PsiH_bar_GHPType{-3, -3};

    GHPType OmH_GHPType{-1, -1};
    GHPType OmH_bar_GHPType{-3, -3};


    HeldScalarsTypeD() = default;

    void fromSpinCoeffs(const SpinCoefficients& sc);
    void print() const;
};

#endif //GHZ_NUMERIC_HELDSCALARS_HPP
