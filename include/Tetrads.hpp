//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_TETRADS_HPP
#define GHZ_NUMERIC_TETRADS_HPP

#pragma once
//#include "Vectors.hpp"
#include "TeukTypes.hpp"
#include "VectorsGHZ.hpp"
#include "SpinCoeffsNP.hpp"
#include "GHPScalars.hpp"
#include "HeldScalars.hpp"
#include "WeylScalars.hpp"
#include "Coords.hpp"
#include "KerrMetric.hpp"

using namespace teuk;

class Tetrad {
protected:
    const KerrMetric& metric;
    const CoordinateSystem& coords;

public:
    ghz::Vector4 l, n;
    ghz::CVector4 m, mbar;
    SpinCoefficients sc;
    SpinCoefficientsGHP sc_ghp;
    HeldCoefficients sc_held;
    WeylScalars weyls;

    explicit Tetrad(KerrMetric& gKerr, CoordinateSystem& c)
            : metric(gKerr), coords(c) {}

    virtual ~Tetrad() = default;
    virtual void build(Real tBL_u_or_v, Real r, Real theta, Real phiBL_in_or_out) = 0;

};

class TetradTransformations {
public:
    static void conformal_transform(Tetrad& tetrad, const Real Om);
    static void spin_boost(Tetrad& tetrad, const Complex& lambda);
    //static void null_rotn_about_l(Tetrad& tetrad, const Complex& c);
    //static void null_rotn_about_n(Tetrad& tetrad, const Complex& c);
};


#endif //GHZ_NUMERIC_TETRADS_HPP
