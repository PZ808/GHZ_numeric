//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_COORDS_HPP
#define GHZ_NUMERIC_COORDS_HPP

#pragma once
#include <cmath>
#include <string>
#include "KerrMetric.hpp"

enum class CoordType {
    BoyerLindquist,
    IngoingKerr,
    OutgoingKerr
};

struct Coords {
    Real x0, x1, x2, x3;

    Coords(Real x0_, Real x1_, Real x2_, Real x3_)
            : x0(x0_), x1(x1_), x2(x2_), x3(x3_) {}
};;


struct BLCoords : Coords {
    using Coords::Coords;
    //double t, r, th, ph;
};

struct IngoingCoords : Coords {
    using Coords::Coords;
    //double v, r, th, ph_in;
};

struct OutgoingCoords : Coords {
    using Coords::Coords;
    //double u, r, z, ph_out; // z = cos(th)
};

class CoordinateSystem {
private:
    const KerrMetric& metric;
    CoordType currentType;

    // helper: analytic tortoise coordinate
    [[nodiscard]] double rStar(double r) const;
    [[nodiscard]] double phiSharp(double r) const;

public:
    explicit CoordinateSystem(const KerrMetric& g, CoordType type = CoordType::BoyerLindquist);

    CoordType type() const { return currentType; }
    void setType(CoordType type);

    // transforms
    IngoingCoords BLtoIngoing(const BLCoords& bl) const;
    OutgoingCoords BLtoOutgoing(const BLCoords& bl) const;
    BLCoords IngoingtoBoyerLindquist(const IngoingCoords& in) const;
    BLCoords OutGoingtoBoyerLindquist(const OutgoingCoords& out) const;

    // print readable form
    std::string typeName() const;
};

#endif //GHZ_NUMERIC_COORDS_HPP
