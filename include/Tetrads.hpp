//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_TETRADS_HPP
#define GHZ_NUMERIC_TETRADS_HPP

#pragma once
#include "Vectors.hpp"
#include "SpinCoeffsNP.hpp"
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

    explicit Tetrad(KerrMetric& gKerr, CoordinateSystem& c)
            : metric(gKerr), coords(c) {}

    virtual ~Tetrad() = default;

    // build the tetrad in the given coordinate system
    virtual void build(double tBL_u_or_v, double r, double theta, double phiBL_in_or_out) = 0;

    // optional: rotate tetrad using a complex null rotation (GHP transformation)
    virtual void rotate(double a_real, double a_imag) {}
    const SpinCoefficients& spinCoeffs() const { return sc; }
};

#endif //GHZ_NUMERIC_TETRADS_HPP
