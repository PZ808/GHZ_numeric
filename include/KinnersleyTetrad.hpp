//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
#define GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
#pragma once

#include "Tetrads.hpp"
// family of Kinnersly tetrad classes — one for every possible CoordT type.
template <typename CoordT>
class KinnersleyTetrad : public Tetrad {
public:

    using Tetrad::Tetrad;
    using Real = teuk::Real;  // or from metric/coords if templated further

    // Coordinate-dependent tetrad builder
    void build_tetrad(const CoordT& X);

    // Common build() signature (dispatches to correct one)
    void build(Real time, Real r, Real polar, Real phi) {
        CoordT X{time, r, polar, phi};
        build_tetrad(X);
    }

    SpinCoefficientsGHP get_spin_coeffs_at(const CoordT& X) const;
    WeylScalars         get_weyl_scalars_at(const CoordT& X) const;
    HeldCoefficients get_held_scalars_at(const CoordT& X) const;
    Scalars get_scalars_at(const CoordT& X) const;
    SpinCoefficientsGHP get_ghp_scalars_at(const OutgoingCoordsCompact &X) const;
};

class KinnersleyTetradBL : public Tetrad {
public:
    using Tetrad::Tetrad;


    void build(Real t_BL, Real r, Real theta, Real ph_BL) override;
    void build_tetrad(const BLCoords& X_bl);
};

class KinnersleyTetradOutgoing : public Tetrad {
public:
    using Tetrad::Tetrad;

    void build(Real u, Real r, Real z, Real ph) override;
    void build_tetrad(const OutgoingCoords& Xout) ;
    void build_tetrad_compact(const OutgoingCoordsCompact& Xout) ;
};


#endif //GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
