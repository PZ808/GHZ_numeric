//
// Created by Peter Zimmerman on 24.10.25.
//

#include "../include/KerrMetric.hpp"
#include "../include/MathMacros.hpp"
#include <cmath>

using std::sqrt;
using std::pow;
using namespace math;

KerrMetric::KerrMetric(const KerrParams& p) : params(p) {
    // set the conformal params
    k2_ = r_minus()/r_plus();
    lambda_ = r_plus();
    mu_ = (1.0+k2_)/2.0;
    alpha_ = k2_;
}

Real KerrMetric::M() const { return params.M; }
Real KerrMetric::a() const { return params.a; }
Real KerrMetric::r_plus() const { return params.M+sqrt(params.M*params.M-params.a*params.a); }
Real KerrMetric::r_minus() const { return params.M-sqrt(params.M*params.M-params.a*params.a); }
Real KerrMetric::Om_plus() const { return params.a/( sqr(params.M+sqrt(params.M*params.M-params.a*params.a))+ sqr(params.a) );  }
Real KerrMetric::Om_minus() const { return params.a / ( sqr(params.M-sqrt(params.M*params.M-params.a*params.a)) + sqr(params.a) );  }
Real KerrMetric::kappa_plus() const { return sqrt(params.M*params.M-params.a*params.a)/(sqr(params.M+sqrt(params.M*params.M - params.a*params.a))+ sqr(params.a) ); }
Real KerrMetric::kappa_minus() const { return sqrt(params.M*params.M-params.a*params.a)/(sqr(params.M-sqrt(params.M*params.M - params.a*params.a))+ sqr(params.a) ) ; }
// conformal quantity getters
Real KerrMetric::k2_C() const {return k2_;}
Real KerrMetric::mu_C() const {return mu_;}
Real KerrMetric::alpha_C() const {return alpha_;}
Real KerrMetric::lambda_C() const {return lambda_;}

Real KerrMetric::Sigma(Real r, Real theta) const {
    Real z = math::Cos(theta);
    return sqr(r) + sqr(z*params.a);
}
Real KerrMetric::Sigma_z(Real r, Real z) const {
    return sqr(r) + sqr(z*params.a);
}
Real KerrMetric::Delta(Real r) const { return sqr(r)- 2.0*params.M*r +sqr(params.a); }

Real KerrMetric::Lambda(Real r, Real theta) const {
    Real z = math::Cos(theta);
    Real s2 = sqr(math::Sin(theta));
    Real del = sqr(r)- two*params.M*r +sqr(params.a);
    return sqr(sqr(r)+sqr(params.a) )-sqr(params.a)*del*s2;
}
Real KerrMetric::Lambda_z(Real r, Real z) const {
    Real s2 =  one - sqr(z);
    Real del = sqr(r)- two*params.M*r +sqr(params.a);
    return sqr(sqr(r)+sqr(params.a) )-sqr(params.a)*del*s2;
}