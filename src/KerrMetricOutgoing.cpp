//
// Created by Peter Zimmerman on 25.10.25.
//

#include <cmath>
#include "../include/KerrMetric.hpp"
#include "../include/KerrMetricOutgoing.hpp"
#include "../include/MathMacros.hpp"
#include <iostream>
using namespace  math;

//
// Xout = u,r,z,phi_out or compacitified u,sigma, z, phi_out
//
KerrMetricOutgoing::KerrMetricOutgoing(const KerrParams& p, const KerrMetric& km) : params(p) , kerr_metric(km) {
    M_ = params.M;
    a_ = params.a;
    lambda_ = kerr_metric.lambda_C();
    rho0_ = kerr_metric.r_plus()/lambda_;
}

// construct local values of the metric in outgoing
// Kerr coordinates
void KerrMetricOutgoing::build(const OutgoingCoords Xout) {
    del_= kerr_metric.Delta(Xout.x1);
    sig_ = sqr(Xout.x1)+sqr(a_*Xout.x2);
    s2_ = 1.0-sqr(Xout.x2);
    s1_ = sqrt(s1_);
}
// compatify
void KerrMetricOutgoing::build_compact_from_outgoing(const OutgoingCoords Xout) {
    Real sigma = lambda_*rho0_/Xout.x1;
    Real r_of_sig = lambda_*rho0_/sigma;
    del_= kerr_metric.Delta(r_of_sig)*sqr(Om_);
    sig_ = sqr(Om_)*kerr_metric.Sigma_z(r_of_sig,Xout.x2);
    s2_ = 1.0-sqr(Xout.x2);
    s1_ = sqrt(s1_);
}
void KerrMetricOutgoing::build_compact(const OutgoingCoordsCompact Xout_C) {
    Real sigma =Xout_C.x1;
    Real r_of_sig = lambda_*rho0_/sigma;
    Om_ = sigma/lambda_;
    del_= kerr_metric.Delta(r_of_sig)*sqr(Om_);
    sig_ = sqr(Om_)*kerr_metric.Sigma_z(r_of_sig,Xout_C.x2);
    s2_ = 1.0-sqr(Xout_C.x2);
    s1_ = sqrt(s1_);
}

ghz::SymmetricMatrix4 KerrMetricOutgoing::g(const OutgoingCoords Xout) const {
    Real r = Xout.x1;

    Real g_uu   = (1.0-2.0*M_*r/sig_);
    Real g_ur   = 1.0;
    Real g_uph  = 2.0*M_*a_*r*s2_/sig_;
    Real g_rph  = -a_*s2_;
    Real g_zz = -sig_/s2_;
    Real g_phph = -( (sqr(r*r+a_*a_)-sqr(a_)*del_*s2_) * s2_)/sig_;

    return { g_uu, g_ur, 0.0, g_uph, 0.0, 0.0, g_rph, g_zz, 0.0, g_phph };
}

ghz::SymmetricMatrix4 KerrMetricOutgoing::g_tilde(const OutgoingCoordsCompact Xout) const {
    Real sigma = Xout.x1;
    Real r_of_sig = lambda_*rho0_/sigma;
    Real r = r_of_sig;

    Real g_uu   = (1.0-2.0*M_*r/sig_) * sqr(Om_);
    Real g_us   = 1.0 * sqr(Om_)*(-lambda_/sqr(sigma));
    Real g_uph  = (2.0*M_*a_*r*s2_/sig_) * sqr(Om_);
    Real g_sph  = -a_*s2_* sqr(Om_)*(-lambda_/sqr(sigma));
    Real g_zz = (-sig_/s2_)* sqr(Om_);
    Real g_phph = ( -( (sqr(r*r+a_*a_)-sqr(a_)*del_*s2_) * s2_ ) / sig_ ) * sqr(Om_);

    return { g_uu, g_us, 0.0, g_uph, 0.0, 0.0, g_sph, g_zz, 0.0, g_phph }  ;
}

ghz::SymmetricMatrix4 KerrMetricOutgoing::ginv(const OutgoingCoords Xout) const {
    Real r = Xout.x1;
    Real z = Xout.x2;


    Real ginv_uu   = 0.0;
    Real ginv_ur   = 0.0;
    Real ginv_uph  = 0.0;
    Real ginv_rr   = 0.0;
    Real ginv_thth = 0.0;
    Real ginv_rph  = 0.0;
    Real ginv_phph = 0.0;

    return { ginv_uu, ginv_ur, 0, ginv_uph, ginv_rr, 0.0, ginv_rph, ginv_thth, 0, ginv_phph };
}