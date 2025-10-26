//
// Created by Peter Zimmerman on 24.10.25.
//

#include "../include/Coords.hpp"

#include <complex>
#include <cmath>
#include <sstream>

using std::sqrt, std::log, std::abs;

CoordinateSystem::CoordinateSystem(const KerrMetric& g, CoordType type)
        : metric(g), currentType(type) {}

void CoordinateSystem::setType(CoordType type) { currentType = type; }

std::string CoordinateSystem::typeName() const {
    switch (currentType) {
        case CoordType::BoyerLindquist: return "BoyerLindquist";
        case CoordType::IngoingKerr:    return "IngoingKerr";
        case CoordType::OutgoingKerr:   return "OutgoingKerr";
        default:                        return "Unknown";
    }
}

// Analytic expression for r_*(r) from Poisson's toolkit pp. 193
double CoordinateSystem::rStar(double r) const {
    double M = metric.M();
    double a = metric.a();
    double r_plus  = metric.r_plus(); // M + sqrt(M*M-a*a);
    double r_minus = metric.r_minus(); //M - sqrt(M*M-a*a);

    //
    double rpm_diff = std::sqrt(M*M-a*a);
    double term1 = ( M*r_plus * log(abs(r/r_plus-1)) ) / rpm_diff;
    double term2 = -( M*r_minus *  log(abs(r/r_minus-1)) ) / rpm_diff;
    return r + term1 + term2;
}

// Analytic expression for phi_sharp(r) from Poisson's toolkit pp. 193
double CoordinateSystem::phiSharp(double r) const {
    double M = metric.M();
    double a = metric.a();
    double r_plus  = metric.r_plus(); // M + sqrt(M*M-a*a);
    double r_minus = metric.r_minus(); //M - sqrt(M*M-a*a);

    // phi_sharp
    double rpm_diff = sqrt(M*M-a*a);
    return (a/(2.*rpm_diff)) * log(abs((r-r_plus)/(r-r_minus)));
}

// ---------------------------
// Transformation definitions
// ---------------------------
IngoingCoords CoordinateSystem::BLtoIngoing(const BLCoords& bl) const {
    // convert from bl to ingoing

    Real t = bl.x0;
    Real r = bl.x1;
    Real th = bl.x2;
    Real ph = bl.x3;

    Real v_in = t + rStar(r);
    Real r_in = r;
    Real th_in = th;
    Real ph_in = ph + phiSharp(r); // optional correction term

    return {v_in, r_in, th_in, ph_in};
}

OutgoingCoords CoordinateSystem::BLtoOutgoing(const BLCoords& bl) const {
    // convert from bl to outgoing
    double u_out  = bl.x0 - rStar(bl.x1);
    double r_out = bl.x1;
    double z_out = std::cos(bl.x2);
    double ph_out = bl.x3 - phiSharp(bl.x1);

    return {u_out, r_out, z_out, ph_out};
}

BLCoords CoordinateSystem::IngoingtoBoyerLindquist(const IngoingCoords &in) const {
    // convert from ingoing to bl
    double t_bl  = in.x0 - rStar(in.x1);
    double r_bl = in.x1;
    double th_bl = in.x2;
    double ph_bl = in.x3 - phiSharp(in.x1);
    return { t_bl, r_bl, th_bl, ph_bl };
}

BLCoords CoordinateSystem::OutGoingtoBoyerLindquist(const OutgoingCoords &out) const {
    // convert from outgoing to bl
    double t_bl  = out.x0 + rStar(out.x1);
    double r_bl = out.x1;
    double th_bl = std::acos(out.x2);
    double ph_bl = out.x3 + phiSharp(r_bl);
    return { t_bl, r_bl, th_bl, ph_bl };
}