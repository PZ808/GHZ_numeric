//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_TETRADS_HPP
#define GHZ_NUMERIC_TETRADS_HPP

#pragma once
#include "Vectors.hpp"
#include "SpinCoeffsNP.hpp"
#include "GHPScalars.hpp"
#include "HeldScalars.hpp"
#include "WeylScalars.hpp"
#include "Coords.hpp"
#include "KerrMetric.hpp"

class Tetrad {
protected:
    const KerrMetric& metric;
    const CoordinateSystem& coords;

public:
    Vector4 l, n;
    CVector4 m, mbar;
    SpinCoefficients sc;
    SpinCoefficientsGHP sc_ghp;
    HeldCoefficients sc_held;
    WeylScalars weyls;

    explicit Tetrad(KerrMetric& gKerr, CoordinateSystem& c)
            : metric(gKerr), coords(c) {}

    virtual ~Tetrad() = default;

    // build the tetrad in the given coordinate system
    virtual void build(double tBL_u_or_v, double r, double theta, double phiBL_in_or_out) = 0;

    //const SpinCoefficients& spinCoeffs() const { return sc; }
    //const GHPScalarsTypeD& ghpScalarsTypeD() const { return ghps; }
    //const SpinCoefficientsGHP& spinCoefficientsGhp() const { return sc_ghp; }
    //const WeylScalars& weylScalars() const { return weyls; }
};

#endif //GHZ_NUMERIC_TETRADS_HPP
