//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_HELDSCALARS_HPP
#define GHZ_NUMERIC_HELDSCALARS_HPP

#pragma once
#include <complex>
#include <string>
#include "GHPScalars.hpp"
#include "WeylScalars.hpp"
#include "SpinCoeffsNP.hpp"
#include "MathMacros.hpp"

using Complex = teuk::Complex;


class HeldScalar : public GHPScalar<Complex> {
    // GHP scalars for type-D spacetime in a Weyl aligned NP frame
    // (shear-free geodesic null congruence)
public:
    void print() const;

};

struct HeldCoefficients {
    GHPScalar<Complex> rhopH, tauH,  rhopH_bar, tauH_bar;
    GHPScalar<Complex> PsiH, PsiH_bar;
    GHPScalar<Complex> OmH, OmH_bar;

    HeldCoefficients() = default;

    HeldCoefficients(const SpinCoefficientsGHP& sc_ghp, const WeylScalars& weyl_scs) {
        // initialize weights according to GHP convention (p,q)
        // (using Held’s sign conventions)
        using SCT = SpinCoeffType;
        Complex rho = sc_ghp.rho.value();
        Complex rhob = std::conj(rho);

        // Held Scalars in Kinnersely tetrad
        rhopH     = GHPScalar<Complex>(-1.0/2.0, -2,-2);
        rhopH_bar = GHPScalar<Complex>(-1.0/2.0, -2, -2);
        tauH      = GHPScalar<Complex>( sc_ghp.tau.value()/(rho*rhob), -1, -3);
        tauH_bar  = GHPScalar<Complex>( std::conj(tauH.value()), -3, -1);
        PsiH      = GHPScalar<Complex>( weyl_scs.get(WeylScalarType::Psi2)/math::cube(rho), -3, -3);
        PsiH_bar  = GHPScalar<Complex>( std::conj(PsiH.value()), -3, -3);
        OmH       = GHPScalar<Complex>((rho-rhob)/(rho*rhob), -1, -1);
        OmH_bar   = GHPScalar<Complex>(std::conj(OmH.value()), -1,-1);
    };

};


#endif //GHZ_NUMERIC_HELDSCALARS_HPP
