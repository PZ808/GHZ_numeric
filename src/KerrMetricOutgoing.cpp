//
// Created by Peter Zimmerman on 25.10.25.
//

#include <cmath>
#include "../include/KerrMetric.hpp"
#include "../include/KerrMetricOutgoing.hpp"
#include "../include/MathMacros.hpp"

using namespace  math;

KerrMetricOutgoing::KerrMetricOutgoing(const KerrParams& p, const KerrMetric& km) : params(p) , kerr_metric(km) {
}

void KerrMetricOutgoing::build(const OutgoingCoords Xout) {
    M_ = params.M;
    a_ = params.a;
    del_= kerr_metric.Delta(Xout.r);
    sig_ = sqr(Xout.r)+sqr(a_*Xout.z);
    s2_ = 1.0-sqr(Xout.z);
    s1_ = sqrt(s1_);
}

ghz::SymmetricMatrix4 KerrMetricOutgoing::g(const OutgoingCoords Xout) const {
    Real r = Xout.r;
    Real z = Xout.z;

    Real M   = params.M;
    Real a   = params.a;
    Real del = kerr_metric.Delta(r);
    Real s2  = 1-sqr(z);
    Real sig = sqr(r)+sqr(a*z); //Sigma(r, th);

    Real g_uu   = (1.0-2.0*M_*r/sig_);
    Real g_ur   = 1.0;
    Real g_uph  = 2.0*M_*a_*r*s2_/sig_;
    Real g_rph  = -a_*s2_;
    Real g_zz = -sig_/s2_;
    Real g_phph = -( (sqr(r*r+a*a)-sqr(a_)*del_*s2_) * s2_)/sig_;

    return { g_uu, g_ur, 0.0, g_uph, 0.0, 0.0, g_rph, g_zz, 0.0, g_phph };
}

ghz::SymmetricMatrix4 KerrMetricOutgoing::ginv(const OutgoingCoords Xout) const {
    Real r = Xout.r;
    Real z = Xout.z;


    Real ginv_uu   = 0.0;
    Real ginv_ur   = 0.0;
    Real ginv_uph  = 0.0;
    Real ginv_rr   = 0.0;
    Real ginv_thth = 0.0;
    Real ginv_rph  = 0.0;
    Real ginv_phph = 0.0;

    return { ginv_uu, ginv_ur, 0, ginv_uph, ginv_rr, 0.0, ginv_rph, ginv_thth, 0, ginv_phph };
}