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
    OutgoingKerr,
    OutgoingKerrCompact,
    Hyperboloidal
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

struct OutgoingCoordsCompact : Coords {
    using Coords::Coords;
    //double u, r, z, ph_out; // z = cos(th)
};

struct HyperboloidalCoords : Coords {
    using Coords::Coords;
    // double tau, sigma, z, phi
};

class CoordinateSystem {
private:
    const KerrMetric& metric;
    CoordType currentType;

    // helper: analytic tortoise coordinate
    [[nodiscard]] Real rStar_(Real r) const;
    [[nodiscard]] Real phiSharp_(Real r) const;
    Real height_(Real sigma) const;
    Real OmegConf_(Real sigma) const;
    Real rho_(Real sigma) const;

public:
    explicit CoordinateSystem(const KerrMetric& g, CoordType type = CoordType::BoyerLindquist);

    CoordType type() const { return currentType; }
    void setType(CoordType type);


    // transforms
    IngoingCoords bl_to_ingoing(const BLCoords& bl) const;
    BLCoords ingoing_to_bl(const IngoingCoords& in) const;
    HyperboloidalCoords ingoing_to_hyperboloidal(const IngoingCoords& in) const;

    [[maybe_unused]] BLCoords outgoing_to_bl(const OutgoingCoords& out) const;

    // print readable form
    std::string typeName() const;

    Real sigma_from_r_(Real r) const;

    Real r_from_sigma(Real r) const;

    Real Omeg_comf(Real sigma) const;

    OutgoingCoords bl_to_outgoing(const BLCoords &bl) const;
};

#endif //GHZ_NUMERIC_COORDS_HPP
