//
// Created by Peter Zimmerman on 25.11.25.
//

#include "ghz/geom/KinnersleyTetrad.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/ghp/GHPScalars.hpp"


HeldCoefficients::HeldCoefficients(const SpinCoefficientsGHP &sc_ghp,
                                                   const WeylScalars &weyl_scs) {
    // initialize weights according to GHP convention (p,q)
    // (using Held’s sign conventions)
    using SCT = SpinCoeffType;
    Complex rho = sc_ghp.rho.value();
    Complex rhob = std::conj(rho);

    // Held Scalars in Kinnersely tetrad
    rhopH     = HeldScalar(teuk::half, -2,-2);
    rhopH_bar = HeldScalar(-teuk::half, -2, -2);
    tauH      = HeldScalar( sc_ghp.tau.value()/(rho*rhob), -1, -3);
    tauH_bar  = HeldScalar( std::conj(tauH.value()), -3, -1);
    PsiH      = HeldScalar( weyl_scs.get(WeylScalarType::Psi2)/math::cube(rho), -3, -3);
    PsiH_bar  = HeldScalar( std::conj(PsiH.value()), -3, -3);
    OmH       = HeldScalar((rho-rhob)/(rho*rhob), -1, -1);
    OmH_bar   = HeldScalar(std::conj(OmH.value()), -1,-1);
};

