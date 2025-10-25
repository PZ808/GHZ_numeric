//
// Created by Peter Zimmerman on 24.10.25.
//

#include "../include/KerrMetric.hpp"
#include "../include/KerrMetricBL.hpp"
#include <cmath>
#include "../include/MathMacros.hpp"

KerrMetricBL::KerrMetricBL(const KerrParams& p, const KerrMetric& km) : params(p) , kerr_metric(km) {
}

//double KerrMetricBL::Sigma(double r, double th) const {
 //   return kerr_metric.Sigma(r,th); //r*r + params.a*params.a*std::cos(th)*std::cos(th);
//}

//double KerrMetricBL::Delta(double r) const {
    //return kerr_metric.Delta(r); // r*r - 2.0*params.M*r + params.a*params.a;
//}

std::array<double, 10> KerrMetricBL::g(double t, double r, double th, double phi) const {
    (void)t; (void)phi;
    // Kerr BL metric in mostly minus
    double s2 = std::sin(th)*std::sin(th);
    double sig =  kerr_metric.Sigma(r,th); //Sigma(r, th);
    double del = kerr_metric.Delta(r);
    double lam = kerr_metric.Lambda(r, th);

    double g_tt     = (1.0 - 2.0*params.M*r/sig);
    double g_rr     = -sig/del;
    double g_thth   = -sig;
    double g_tphi   = 2.0*params.M*params.a*r*s2/sig;
    double g_phiphi = -lam*s2/sig;

    return { g_tt, 0, 0, g_tphi, g_rr, 0, 0, g_thth, 0, g_phiphi };
}

std::array<double, 10> KerrMetricBL::ginv(double t, double r, double th, double phi) const {
    (void)t; (void)phi;
    double s2 = std::sin(th)*std::sin(th);
    double sig = kerr_metric.Sigma(r, th);
    double del = kerr_metric.Delta(r);
    double lam = kerr_metric.Lambda(r, th);

    double ginv_tt     = -lam/(sig*del);
    double ginv_rr     = del/sig;
    double ginv_thth   = 1/sig;
    double ginv_tphi   = -2.*params.M*params.a*r/(sig*del);
    double ginv_phiphi = (del-math::sqr(params.a)*s2)/(sig*del*s2);

    return { ginv_tt, 0, 0, ginv_tphi, ginv_rr, 0, 0, ginv_thth, 0, ginv_phiphi };
}