//
// Created by Peter Zimmerman on 24.10.25.
//

#include "../include/Coords.hpp"

#include <complex>
#include <cmath>
#include <sstream>
#include "../include/MathMacros.hpp"

using std::sqrt, std::log, std::abs;
using math::sqr;
using namespace teuk;

CoordinateSystem::CoordinateSystem(const KerrMetric& g, CoordType type)
        : metric(g), currentType(type) {}

void CoordinateSystem::setType(CoordType type) { currentType = type; }

std::string CoordinateSystem::typeName() const {
    switch (currentType) {
        case CoordType::BoyerLindquist: return "BoyerLindquist";
        case CoordType::IngoingKerr:    return "IngoingKerr";
        case CoordType::OutgoingKerr:   return "OutgoingKerr";
        case CoordType::OutgoingKerrCompact:   return "OutgoingKerr";
        case CoordType::Hyperboloidal:  return "Hyperboloidal";
        default:                        return "Unknown";
    }
}

// Analytic expression for r_*(r) from Poisson's toolkit pp. 193
Real CoordinateSystem::rStar_(Real r) const {
    Real M = metric.M();
    Real a = metric.a();
    Real r_plus  = metric.r_plus(); // M + sqrt(M*M-a*a);
    Real r_minus = metric.r_minus(); //M - sqrt(M*M-a*a);

    //
    Real rpm_diff = std::sqrt(M*M-a*a);
    Real term1 = ( M*r_plus * log(abs(r/r_plus-1)) ) / rpm_diff;
    Real term2 = -( M*r_minus *  log(abs(r/r_minus-1)) ) / rpm_diff;
    return r + term1 + term2;
}

// Analytic expression for phi_sharp(r) from Poisson's toolkit pp. 193
Real CoordinateSystem::phiSharp_(Real r) const {
    Real M = metric.M();
    Real a = metric.a();
    Real r_plus  = metric.r_plus(); // M + sqrt(M*M-a*a);
    Real r_minus = metric.r_minus(); //M - sqrt(M*M-a*a);

    // phi_sharp
    Real rpm_diff = sqrt(M*M-a*a);
    return (a/(2.*rpm_diff)) * log(abs((r-r_plus)/(r-r_minus)));
}

Real CoordinateSystem::sigma_from_r_(Real r) const {
    return metric.lambda_C()/r;
}
Real CoordinateSystem::r_from_sigma(Real r) const {
    return r/metric.lambda_C();
}

Real CoordinateSystem::height_(Real sigma) const {
    // minimal gauge
    Real M = metric.M();
    Real lambda = metric.lambda_C();
    Real k2 = metric.k2_C();
    Real mu = metric.mu_C();
    Real rho0 = 1.0; // enforced by minimal gauge condition and choice of k2,mu,and lambda
    return -2.0*rho0*(1.0/sigma-2.0*mu*std::log(sigma)/rho0);
}

Real CoordinateSystem::Omeg_comf(teuk::Real sigma) const {
    Real lambda = metric.r_plus();
    return sigma/lambda;
}
// ---------------------------
// Transformation definitions
// ---------------------------
IngoingCoords CoordinateSystem::bl_to_ingoing(const BLCoords& bl) const {
    // convert from bl to ingoing

    Real t = bl.x0;
    Real r = bl.x1;
    Real th = bl.x2;
    Real ph = bl.x3;

    Real v_in = t + rStar_(r);
    Real r_in = r;
    Real th_in = th;
    Real ph_in = ph + phiSharp_(r); // optional correction term

    return {v_in, r_in, th_in, ph_in};
}

OutgoingCoords CoordinateSystem::bl_to_outgoing(const BLCoords& bl) const {
    // convert from bl to outgoing
    double u_out  = bl.x0 - rStar_(bl.x1);
    double r_out = bl.x1;
    double z_out = std::cos(bl.x2);
    double ph_out = bl.x3 - phiSharp_(bl.x1);

    return {u_out, r_out, z_out, ph_out};
}

BLCoords CoordinateSystem::ingoing_to_bl(const IngoingCoords &in) const {
    // convert from ingoing to bl
    double t_bl  = in.x0 - rStar_(in.x1);
    double r_bl = in.x1;
    double th_bl = in.x2;
    double ph_bl = in.x3 - phiSharp_(in.x1);
    return { t_bl, r_bl, th_bl, ph_bl };
}

[[maybe_unused]] BLCoords CoordinateSystem::outgoing_to_bl(const OutgoingCoords &out) const {
    // convert from outgoing to bl
    double t_bl  = out.x0 + rStar_(out.x1);
    double r_bl = out.x1;
    double th_bl = std::acos(out.x2);
    double ph_bl = out.x3 + phiSharp_(r_bl);
    return { t_bl, r_bl, th_bl, ph_bl };
}