//
// Created by Peter Zimmerman on 24.10.25.
//

#include "ghz/geom/KerrMetric.hpp"
#include "ghz/core/MathMacros.hpp"
#include "ghz/geom/KerrMetricBL.hpp"
#include "ghz/core/MathMacros.hpp"

#include <cassert>
#include <cmath>

using namespace  math;
Complex I = teuk::I;


KerrMetricBL::KerrMetricBL(const KerrParams& p, const KerrMetric& km)
        : params(p) , kerr_metric(km) {
        // No cached values yet
        cache_valid_ = false;
}
void KerrMetricBL::build_at(const BLCoords Xbl) {


    // Cache the coordinate point and mark valid
    X_cached_ = Xbl;
    cache_valid_ = true;
    // Compute all cached scalars at point Xbl
    M_ = params.M;
    a_ = params.a;
    del_= kerr_metric.Delta(Xbl.x1);
    sig_ = sqr(Xbl.x1)+sqr(a_*Xbl.x2);
    s2_ = teuk::one-sqr(Xbl.x2);
    s1_ = sqrt(s1_);
}

teuk::SymmetricMatrix4 KerrMetricBL::g(const BLCoords Xbl) const {

    assert(cache_valid_ && "KerrMetricBL::g() called before build()");
    // Kerr BL metric in mostly minus
    Real r = Xbl.x1; // r = r
    Real z = Xbl.x2; // z = \cos\theta
    Real lam = kerr_metric.Lambda_z(r, z);

    // Construct matrix using cached sig_, del_, etc.
    Real g_tt     = (one - two*M_*r/sig_);
    Real g_rr     = -sig_/del_;
    Real g_zz   = -sig_/s2_;
    Real g_tphi   = two * M_*a_*r*s2_/sig_;
    Real g_phiphi = -lam*s2_/sig_;

    return { g_tt, 0, 0, g_tphi, g_rr, 0, 0, g_zz, 0, g_phiphi };
}

teuk::SymmetricMatrix4 KerrMetricBL::ginv(const BLCoords Xbl) const {

    assert(cache_valid_ && "KerrMetricBL::g() called before build()");

    Real r = Xbl.x1;
    Real th = Xbl.x2;
    Real lam = kerr_metric.Lambda(r, th);

    Real ginv_tt     = -lam/(sig_*del_);
    Real ginv_rr     = del_/sig_;
    Real ginv_zz   = s2_/sig_;
    Real ginv_tphi   = -two*params.M*params.a*r/(sig_*del_);
    Real ginv_phiphi = (del_-math::sqr(params.a)*s2_)/(sig_*del_*s2_);

    return { ginv_tt, 0, 0, ginv_tphi, ginv_rr, 0, 0, ginv_zz, 0, ginv_phiphi };
}
