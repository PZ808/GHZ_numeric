//
// Created by Peter Zimmerman on 25.11.25.
//

#include "ghz/geom/KinnersleyTetrad.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/ghp/GHPScalars.hpp"

#include <algorithm>


ghp::HeldCoefficients::HeldCoefficients(const teuk::Real mass,
                                        const teuk::Real spin,
                                        const teuk::Real z) {
    const teuk::Real sin_theta = math::Sqrt(
            std::max(teuk::Real(0), teuk::Real(1) - z*z));

    // Analytic Held scalars for the Kerr Kinnersley tetrad with signature
    // (+---).  Keeping these expressions independent of rho avoids forming
    // small radial quantities and dividing them back out.
    rhopH     = HeldScalar(-teuk::half, -2, -2);
    rhopH_bar = HeldScalar(-teuk::half, -2, -2);
    tauH      = HeldScalar(-teuk::I * spin * sin_theta /
                           math::Sqrt(teuk::Real(2)), -1, -3);
    tauH_bar  = HeldScalar(std::conj(tauH.value()), -3, -1);
    PsiH      = HeldScalar(Complex(mass, teuk::Real(0)), -3, -3);
    PsiH_bar  = HeldScalar(Complex(mass, teuk::Real(0)), -3, -3);
    OmH       = HeldScalar(-teuk::Real(2) * teuk::I * spin * z, -1, -1);
    OmH_bar   = HeldScalar(std::conj(OmH.value()), -1, -1);
}


ghp::HeldCoefficients::HeldCoefficients(const ghp::SpinCoefficientsGHP &sc_ghp,
                                                   const WeylScalars &weyl_scs) {
    // initialize weights according to GHP convention (p,q)
    // (using Held’s sign conventions)
    Complex rho = sc_ghp.rho.value();
    Complex rhob = std::conj(rho);

    // Generic definitions.  Kerr callers should prefer the analytic
    // (mass, spin, z) constructor above.
    rhopH     = HeldScalar(-teuk::half, -2, -2);
    rhopH_bar = HeldScalar(-teuk::half, -2, -2);
    tauH      = HeldScalar( sc_ghp.tau.value()/(rho*rhob), -1, -3);
    tauH_bar  = HeldScalar( std::conj(tauH.value()), -3, -1);
    PsiH      = HeldScalar( weyl_scs.get(WeylScalarType::Psi2)/math::cube(rho), -3, -3);
    PsiH_bar  = HeldScalar( std::conj(PsiH.value()), -3, -3);
    OmH       = HeldScalar((rho-rhob)/(rho*rhob), -1, -1);
    OmH_bar   = HeldScalar(std::conj(OmH.value()), -1,-1);
};
