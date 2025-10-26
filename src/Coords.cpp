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
IngoingCoords CoordinateSystem::toIngoing(const BLCoords& bl) const {
    // convert from bl to ingoing

    double v_in  = bl.t + rStar(bl.r);
    double r_in = bl.r;
    double th_in = bl.th;
    double ph_in = bl.ph + phiSharp(bl.r); // optional correction term

    return {v_in, r_in, th_in, ph_in};
}

OutgoingCoords CoordinateSystem::toOutgoing(const BLCoords& bl) const {
    // convert from bl to outgoing
    double u_out  = bl.t - rStar(bl.r);
    double r_out = bl.r;
    double z_out = std::cos(bl.th);
    double ph_out = bl.ph - phiSharp(bl.r);

    return {u_out, r_out, z_out, ph_out};
}

BLCoords CoordinateSystem::toBoyerLindquist(const IngoingCoords& in) const {
    // convert from ingoing to bl
    double t_bl  = in.v - rStar(in.r);
    double r_bl = in.r;
    double th_bl = in.th;
    double ph_bl = in.ph_in - phiSharp(r_bl);
    return { t_bl, r_bl, th_bl, ph_bl };
}

BLCoords CoordinateSystem::toBoyerLindquist(const OutgoingCoords& out) const {
    // convert from outgoing to bl
    double t_bl  = out.u + rStar(out.r);
    double r_bl = out.r;
    double th_bl = std::acos(out.z);
    double ph_bl = out.ph_out + phiSharp(r_bl);
    return { t_bl, r_bl, th_bl, ph_bl };
}