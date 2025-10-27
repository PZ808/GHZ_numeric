//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
#define GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
#pragma once

#include "Tetrads.hpp"

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
