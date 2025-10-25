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

struct BLCoords {
    double t, r, th, ph;
};

struct IngoingCoords {
    double v, r, th, ph_in;
};

struct OutgoingCoords {
    double u, r, th, ph_out;
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
    IngoingCoords toIngoing(const BLCoords& bl) const;
    OutgoingCoords toOutgoing(const BLCoords& bl) const;
    BLCoords toBoyerLindquist(const IngoingCoords& in) const;
    BLCoords toBoyerLindquist(const OutgoingCoords& out) const;

    // print readable form
    std::string typeName() const;
};

#endif //GHZ_NUMERIC_COORDS_HPP
