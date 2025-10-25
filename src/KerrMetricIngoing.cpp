//
// Created by Peter Zimmerman on 25.10.25.
//


#include "../include/KerrMetric.hpp"
#include "../include/KerrMetricIngoing.hpp"
#include <cmath>
#include "../include/MathMacros.hpp"

KerrMetricIngoing::KerrMetricIngoing(const KerrParams& p, const KerrMetric& km) : params(p) , kerr_metric(km) {
}

std::array<double, 10> KerrMetricIngoing::g(double v, double r, double th, double phi_in) const {
    (void)v; (void)phi_in;
    double M = params.M;
    double a = params.a;
    double s2 = std::sin(th)*std::sin(th);
    double sig =  kerr_metric.Sigma(r,th); //Sigma(r, th);
    double del = kerr_metric.Delta(r);
    double lam = kerr_metric.Lambda(r, th);

    double g_vv   = -(1.-2.*M*r/sig);
    double g_vr   = 1.0;
    double g_vph   = -a*s2*(1 - 2*M*r/sig);
    double g_rph   = -a*s2;
    double g_thth   = sig;
    double g_phph   = (r*r + a*a + 2*M*r*a*a*s2/sig)*s2;

    return { g_vv, g_vr, 0.0, g_vph, 0.0, 0.0, g_rph, g_thth, 0.0, g_phph };
}

std::array<double, 10> KerrMetricIngoing::ginv(double t, double r, double th, double phi) const {
    (void)t; (void)phi;
    double s2 = std::sin(th)*std::sin(th);
    double sig = kerr_metric.Sigma(r, th);
    double del = kerr_metric.Delta(r);
    double lam = kerr_metric.Lambda(r, th);

    double ginv_vv    = 0.0;
    double ginv_vr    = 0.0;
    double ginv_rr    = 0.0;
    double ginv_rph   = 0.0;
    double ginv_thth  = 0.0;
    double ginv_vph   = 0.0;
    double ginv_phiph = 0.0;

    return { ginv_vv, ginv_vr, 0, ginv_vph, ginv_rr, 0, ginv_rph, ginv_thth, 0, ginv_phiph };
}