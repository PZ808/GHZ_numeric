//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_TETRADS_HPP
#define GHZ_NUMERIC_TETRADS_HPP

#pragma once
//#include "Vectors.hpp"
//#include "TeukTypes.hpp"
#include "ghz/core/GhzTypes.hpp"
#include "ghz/core/VectorsGHZ.hpp"
#include "ghz/ghp/SpinCoeffsNP.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/ghp/WeylScalars.hpp"
#include "ghz/geom/Coords.hpp"
#include "ghz/geom/KerrMetric.hpp"

using namespace teuk;

class Tetrad {
protected:
    const KerrMetric& metric;
    const CoordinateHelper& coord_helper;

public:
    teuk::Vector4 l, n;
    teuk::CVector4 m, mbar;
    SpinCoefficients sc;
    ghp::SpinCoefficientsGHP sc_ghp;
    ghp::HeldCoefficients sc_held;
    ghp::WeylScalars weyls;

    struct Scalars {
        ghp::SpinCoefficientsGHP ghp_scalars;
        ghp::HeldCoefficients held_scalars;
    };

    explicit Tetrad(KerrMetric& gKerr, CoordinateHelper& ch)
            : metric(gKerr), coord_helper(ch) {}

    virtual ~Tetrad() = default;

};

class TetradTransformations {
public:
    static void conformal_transform(Tetrad& tetrad, const Real Om);
    static void spin_boost(Tetrad& tetrad, const Complex& lambda);
    //static void null_rotn_about_l(Tetrad& tetrad, const Complex& c);
    //static void null_rotn_about_n(Tetrad& tetrad, const Complex& c);
};


#endif //GHZ_NUMERIC_TETRADS_HPP
