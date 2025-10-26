//
// Created by Peter Zimmerman on 24.10.25.
//

#include "../include/KerrMetric.hpp"
#include "../include/MathMacros.hpp"
#include "../include/KerrMetricBL.hpp"
#include <cmath>
#include "../include/MathMacros.hpp"

using namespace  math;

KerrMetricBL::KerrMetricBL(const KerrParams& p, const KerrMetric& km) : params(p) , kerr_metric(km) {
}
void KerrMetricBL::build(const BLCoords Xbl) {
    M_ = params.M;
    a_ = params.a;
    del_= kerr_metric.Delta(Xbl.x1);
    sig_ = sqr(Xbl.x1)+sqr(a_*Xbl.x2);
    s2_ = 1.0-sqr(Xbl.x2);
    s1_ = sqrt(s1_);
}

ghz::SymmetricMatrix4 KerrMetricBL::g(const BLCoords Xbl) const {
    // Kerr BL metric in mostly minus
    Real r = Xbl.x1;
    Real th = Xbl.x2;
    double lam = kerr_metric.Lambda(r, th);

    double g_tt     = (1.0 - 2.0*M_*r/sig_);
    double g_rr     = -sig_/del_;
    double g_thth   = -sig_;
    double g_tphi   = 2.0 * M_*a_*r*s2_/sig_;
    double g_phiphi = -lam*s2_/sig_;

    return { g_tt, 0, 0, g_tphi, g_rr, 0, 0, g_thth, 0, g_phiphi };
}

ghz::SymmetricMatrix4 KerrMetricBL::ginv(const BLCoords Xbl) const {
    Real r = Xbl.x1;
    Real th = Xbl.x2;
    double lam = kerr_metric.Lambda(r, th);

    double ginv_tt     = -lam/(sig_*del_);
    double ginv_rr     = del_/sig_;
    double ginv_thth   = 1/sig_;
    double ginv_tphi   = -2.*params.M*params.a*r/(sig_*del_);
    double ginv_phiphi = (del_-math::sqr(params.a)*s2_)/(sig_*del_*s2_);

    return { ginv_tt, 0, 0, ginv_tphi, ginv_rr, 0, 0, ginv_thth, 0, ginv_phiphi };
}