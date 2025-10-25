//
// Created by Peter Zimmerman on 25.10.25.
//

#include <cmath>
#include "../include/KerrMetric.hpp"
#include "../include/KerrMetricOutgoing.hpp"
#include "../include/MathMacros.hpp"

KerrMetricOutgoing::KerrMetricOutgoing(const KerrParams& p, const KerrMetric& km) : params(p) , kerr_metric(km) {
}

std::array<double, 10> KerrMetricOutgoing::g(double u, double r, double th, double phi_out) const {
    (void)u; (void)phi_out;
    double M   = params.M;
    double a   = params.a;
    double s2  = std::sin(th)*std::sin(th);
    double sig =  kerr_metric.Sigma(r,th); //Sigma(r, th);
    double del = kerr_metric.Delta(r);
    double lam = kerr_metric.Lambda(r, th);

    double g_uu   = -(1.-2.*M*r/sig);
    double g_ur   = -1.0;
    double g_uph  = a*s2*(1 - 2*M*r/sig);
    double g_rph  = -a*s2;
    double g_thth = sig;
    double g_phph = (r*r + a*a + 2*M*r*a*a*s2/sig)*s2;

    return { g_uu, g_ur, 0.0, g_uph, 0.0, 0.0, g_rph, g_thth, 0.0, g_phph };
}

std::array<double, 10> KerrMetricOutgoing::ginv(double t, double r, double th, double phi) const {
    (void)t; (void)phi;
    double s2  = std::sin(th)*std::sin(th);
    double sig = kerr_metric.Sigma(r, th);
    double del = kerr_metric.Delta(r);
    double lam = kerr_metric.Lambda(r, th);

    double ginv_uu   = 0.0;
    double ginv_ur   = 0.0;
    double ginv_uph  = 0.0;
    double ginv_rr   = 0.0;
    double ginv_thth = 0.0;
    double ginv_rph  = 0.0;
    double ginv_phph = 0.0;

    return { ginv_uu, ginv_ur, 0, ginv_uph, ginv_rr, 0.0, ginv_rph, ginv_thth, 0, ginv_phph };
}