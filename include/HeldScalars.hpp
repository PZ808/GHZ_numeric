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
#include "SpectralGHPField.hpp"

using Complex = teuk::Complex;


class HeldScalar : public GHPScalar<Complex> {
public:
    HeldScalar() = default;
    HeldScalar(const Complex& v, int p, int q) : GHPScalar<Complex>(v, p, q) {}

};

class HeldField : public GHPField {
public:
    HeldField() = default;
    HeldField(int Nz,
              Scalar init = Scalar(teuk::zeroC, 0, 0), int p = 0, int q = 0)
:     GHPField(1, Nz, init, p, q) {}

    Scalar& operator()(int iz) {
        return GHPField::operator()(0, iz);
    }
    const Scalar& operator()(int iz) const {
        return GHPField::operator()(0, iz);
    }

};

struct HeldCoefficients {
    HeldCoefficients(const SpinCoefficientsGHP &sc_ghp, const WeylScalars &weyl_scs);

    HeldScalar rhopH, tauH,  rhopH_bar, tauH_bar;
    HeldScalar PsiH, PsiH_bar;
    HeldScalar OmH, OmH_bar;

    HeldCoefficients() = default;

};

struct HeldBackgroundFields {
    HeldField rhopH, tauH, tauH_bar, PsiH, PsiH_bar, OmH, OmH_bar;
    HeldBackgroundFields(int Nz) :
            rhopH(Nz, HeldScalar(teuk::zeroC,-2,-2), -2,-2),
            tauH (Nz, HeldScalar(teuk::zeroC,-1,-3), -1,-3),
            tauH_bar(Nz, HeldScalar(teuk::zeroC,-3,-1), -3,-1),
            PsiH (Nz, HeldScalar(teuk::zeroC,-3,-3), -3,-3),
            PsiH_bar(Nz, HeldScalar(teuk::zeroC,-3,-3), -3,-3),
            OmH  (Nz, HeldScalar(teuk::zeroC,-1,-1), -1,-1),
            OmH_bar(Nz,HeldScalar(teuk::zeroC,-1,-1), -1,-1) {}
};

template <typename TetradType, typename CoordType>
HeldBackgroundFields build_held_fields(TetradType& tetrad,
                             const std::vector<teuk::Real>& z_nodes,
                             CoordType &X)
{
    const int Nz = static_cast<int>(z_nodes.size());
    HeldBackgroundFields held_fields(Nz);

    for (int iz = 0; iz < Nz; ++iz) {
        X.x2 = z_nodes[iz];       // set current z to LGL collocation point

        tetrad.build_tetrad(X);   // build tetrad at this point

        auto scalars = tetrad.get_scalars_at(X);  // get all scalars
        // assign directly to HeldFields
        held_fields.rhopH(iz)     = scalars.held_scalars.rhopH;
        held_fields.tauH(iz)      = scalars.held_scalars.tauH;
        held_fields.tauH_bar(iz)  = scalars.held_scalars.tauH_bar;
        held_fields.PsiH(iz)      = scalars.held_scalars.PsiH;
        held_fields.PsiH_bar(iz)  = scalars.held_scalars.PsiH_bar;
        held_fields.OmH(iz)       = scalars.held_scalars.OmH;
        held_fields.OmH_bar(iz)   = scalars.held_scalars.OmH_bar;
    }

    return held_fields;
}

#endif //GHZ_NUMERIC_HELDSCALARS_HPP
