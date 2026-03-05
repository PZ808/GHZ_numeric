//
// Created by Peter Zimmerman on 24.10.25.
//

#include "ghz/geom/Coords.hpp"
#include "ghz/core/MathMacros.hpp"

#include <complex>
#include <cmath>
#include <sstream>

using math::Sqrt, math::Log, math::Abs;
using namespace teuk;

CoordinateHelper::CoordinateHelper(const KerrMetric& g, CoordType type)
        : metric(g), currentType(type) {}

void CoordinateHelper::setType(CoordType type) { currentType = type; }

std::string CoordinateHelper::typeName() const {
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
Real CoordinateHelper::rStar_(Real r) const {
    Real M = metric.M();
    Real a = metric.a();
    Real r_plus  = metric.r_plus(); // M + sqrt(M*M-a*a);
    Real r_minus = metric.r_minus(); //M - sqrt(M*M-a*a);

    //
    Real rpm_diff = math::Sqrt(M*M-a*a);
    Real term1 = ( M*r_plus * log(abs(r/r_plus-1)) ) / rpm_diff;
    Real term2 = -( M*r_minus *  log(abs(r/r_minus-1)) ) / rpm_diff;
    return r + term1 + term2;
}

// Analytic expression for phi_sharp(r) from Poisson's toolkit pp. 193
Real CoordinateHelper::phiSharp_(Real r) const {
    Real M = metric.M();
    Real a = metric.a();
    Real r_plus  = metric.r_plus(); // M + sqrt(M*M-a*a);
    Real r_minus = metric.r_minus(); //M - sqrt(M*M-a*a);

    // phi_sharp
    Real rpm_diff = sqrt(M*M-a*a);
    return (a/(2.*rpm_diff)) * log(abs((r-r_plus)/(r-r_minus)));
}

Real CoordinateHelper::sigma_from_r_(Real r) const {
    return metric.lambda_C()/r;
}
Real CoordinateHelper::r_from_sigma(Real r) const {
    return r/metric.lambda_C();
}

Real CoordinateHelper::height_(Real sigma) const {
    // minimal gauge
    Real M = metric.M();
    Real lambda = metric.lambda_C();
    Real k2 = metric.k2_C();
    Real mu = metric.mu_C();
    Real rho0 = 1.0; // enforced by minimal gauge condition and choice of k2,mu,and lambda
    return -2.0*rho0*(1.0/sigma-2.0*mu*Log(sigma)/rho0);
}

Real CoordinateHelper::Omeg_comf(teuk::Real sigma) const {
    Real lambda = metric.r_plus();
    return sigma/lambda;
}
// ---------------------------
// Transformation definitions
// ---------------------------
IngoingCoords CoordinateHelper::bl_to_ingoing(const BLCoords& bl) const {
    // convert from bl to ingoing

    Real t = bl.x0;
    Real r = bl.x1;
    Real z = bl.x2;
    Real ph = bl.x3;

    Real v_in = t + rStar_(r);
    Real r_in = r;
    Real z_in = z;
    Real ph_in = ph + phiSharp_(r); // optional correction term

    return {v_in, r_in, z_in, ph_in};
}

OutgoingCoords CoordinateHelper::bl_to_outgoing(const BLCoords& bl) const {
    // convert from bl to outgoing
    Real u_out  = bl.x0 - rStar_(bl.x1);
    Real r_out = bl.x1;
    Real z_out = math::Cos(bl.x2);
    Real ph_out = bl.x3 - phiSharp_(bl.x1);

    return {u_out, r_out, z_out, ph_out};
}

BLCoords CoordinateHelper::ingoing_to_bl(const IngoingCoords &in) const {
    // convert from ingoing to bl
    Real t_bl  = in.x0 - rStar_(in.x1);
    Real r_bl = in.x1;
    Real z_bl = in.x2;
    Real ph_bl = in.x3 - phiSharp_(in.x1);
    return { t_bl, r_bl, z_bl, ph_bl };
}

[[maybe_unused]] BLCoords CoordinateHelper::outgoing_to_bl(const OutgoingCoords &out) const {
    // convert from outgoing to bl
    Real t_bl  = out.x0 + rStar_(out.x1);
    Real r_bl = out.x1;
    Real z_bl = acos(out.x2);
    Real ph_bl = out.x3 + phiSharp_(r_bl);
    return { t_bl, r_bl, z_bl, ph_bl };
}