//
// Created by Peter Zimmerman on 25.10.25.
//


#include "../include/KerrMetric.hpp"
#include "../include/KerrMetricIngoing.hpp"
#include <cmath>
#include "../include/MathMacros.hpp"
using namespace  math;

KerrMetricIngoing::KerrMetricIngoing(const KerrParams& p, const KerrMetric& km) : params(p) , kerr_metric(km) {
}
void KerrMetricIngoing::build(const OutgoingCoords Xout) {
    M_ = params.M;
    a_ = params.a;
    del_= kerr_metric.Delta(Xout.x1);
    sig_ = sqr(Xout.x1)+sqr(a_*Xout.x2);
    s2_ = 1.0-sqr(Xout.x2);
    s1_ = sqrt(s1_);
}

ghz::SymmetricMatrix4 KerrMetricIngoing::g(const IngoingCoords Xin) const {
    Real r = Xin.x1;

    Real g_vv   = (1.-2.*M_*r/sig_);
    Real g_vr   = -1.0;
    Real g_vph   = a_*s2_*(1.0-2.0*M_*r/sig_);
    Real g_rph   = a_*s2_;
    Real g_thth   = -sig_;
    Real g_phph   = -(r*r + a_*a_ + 2.0*M_*r*a_*a_*s2_/sig_)*s2_;

    return { g_vv, g_vr, 0.0, g_vph, 0.0, 0.0, g_rph, g_thth, 0.0, g_phph };
}

ghz::SymmetricMatrix4 KerrMetricIngoing::ginv(const IngoingCoords Xin) const {

    Real ginv_vv    = 0.0;
    Real ginv_vr    = 0.0;
    Real ginv_rr    = 0.0;
    Real ginv_rph   = 0.0;
    Real ginv_thth  = 0.0;
    Real ginv_vph   = 0.0;
    Real ginv_phiph = 0.0;

    return { ginv_vv, ginv_vr, 0, ginv_vph, ginv_rr, 0, ginv_rph, ginv_thth, 0, ginv_phiph };
}