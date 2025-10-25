//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
#define GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
#pragma once
#include "SpinCoeffsNP.hpp"
#include "Tetrads.hpp"

class KinnersleyTetradBL : public Tetrad {
public:
    using Tetrad::Tetrad;

    void build(double t_BL, double r, double theta, double ph_BL) override;
};

class KinnersleyTetradOutgoing : public Tetrad {
public:
using Tetrad::Tetrad;

void build(double u, double r, double theta, double ph_out) override;
};

#endif //GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
