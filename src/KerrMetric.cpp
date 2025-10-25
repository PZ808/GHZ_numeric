//
// Created by Peter Zimmerman on 24.10.25.
//

#include "../include/KerrMetric.hpp"
#include "../include/MathMacros.hpp"
#include <cmath>

using std::sqrt;
using std::pow;
using namespace math;

KerrMetric::KerrMetric(const KerrParams& p) : params(p) {}

double KerrMetric::M() const { return params.M; }
double KerrMetric::a() const { return params.a; }
double KerrMetric::r_plus() const { return params.M+sqrt(params.M*params.M-params.a*params.a); }
double KerrMetric::r_minus() const { return params.M-sqrt(params.M*params.M-params.a*params.a); }
double KerrMetric::Om_plus() const { return params.a/( sqr(params.M+sqrt(params.M*params.M-params.a*params.a))+ sqr(params.a) );  }
double KerrMetric::Om_minus() const { return params.a / ( sqr(params.M-sqrt(params.M*params.M-params.a*params.a)) + sqr(params.a) );  }

double KerrMetric::kappa_plus() const { return sqrt(params.M*params.M-params.a*params.a)/(sqr(params.M+sqrt(params.M*params.M - params.a*params.a))+ sqr(params.a) ); }
double KerrMetric::kappa_minus() const { return sqrt(params.M*params.M-params.a*params.a)/(sqr(params.M-sqrt(params.M*params.M - params.a*params.a))+ sqr(params.a) ) ; }


double KerrMetric::Sigma(double r, double theta) const {
    double z = std::cos(theta);
    return sqr(r) + sqr(z*params.a);
}
double KerrMetric::Delta(double r) const {
    return sqr(r)- 2.0*params.M*r +sqr(params.a);
}
double KerrMetric::Lambda(double r, double theta) const {
    double z = std::cos(theta);
    double s2 = sqr(std::sin(theta));
    double del = sqr(r)- 2.0*params.M*r +sqr(params.a);
    return sqr(sqr(r)+sqr(params.a) )-sqr(params.a)*del*s2;
}


#if 0
std::array<double,10> KerrMetric::g(double t, double r, double th, double ph) const {
    (void)t; (void)ph;


    double s2 = std::sin(th)*std::sin(th);
    double sig = Sigma(r, th);
    double del = Delta(r);


    double g_tt = -(1.0 - 2.0*params.M*r/sig);
    double g_rr = sig/del;
    double g_thth = sig;
    double g_phph = (r*r + params.a*params.a + 2.0*params.M*r*params.a*params.a*s2/sig) * s2;
    double g_tph = -2.0*params.M*r*params.a*s2/sig;


    return { g_tt, 0.0, 0.0, g_tph, g_rr, 0.0, 0.0, g_thth, 0.0, g_phph };
}
#endif