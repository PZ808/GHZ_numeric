//
// Created by Peter Zimmerman on 24.10.25.
//

#include "../include/KerrMetric.hpp"
#include "../include/MathMacros.hpp"
#include <cmath>

using std::sqrt;
using std::pow;
using namespace math;

/*
 *
 */
KerrMetric::KerrMetric(const KerrParams& p) : params(p) {
    //
    // set the conformal params
    // (see https://arxiv.org/pdf/1910.13452)
    k2_ = r_minus()/r_plus(); // \kappa^2 := r_-/r_+ , \kappa := j/(1+\sqrt{1-j^2}) j = a/M
    lambda_ = r_plus();  // preferred length scale
    // dimensionless mass and spin parameters
    mu_ = (1.0+k2_)/2.0; // \mu := M /\lambda
    alpha_ = k2_;        // \alpha := a/\lambda
}
//
// Kerr param getters
//
Real KerrMetric::M() const { return params.M; } // mass
Real KerrMetric::a() const { return params.a; } // spin
Real KerrMetric::r_plus() const { return params.M+sqrt(params.M*params.M-params.a*params.a); } // outter horizon
Real KerrMetric::r_minus() const { return params.M-sqrt(params.M*params.M-params.a*params.a); } // inner horizon
Real KerrMetric::Om_plus() const { return params.a/(sqr(params.M+sqrt(params.M*params.M-params.a*params.a))+ sqr(params.a) );  }
Real KerrMetric::Om_minus() const { return params.a / ( sqr(params.M-sqrt(params.M*params.M-params.a*params.a)) + sqr(params.a) );  }
Real KerrMetric::kappa_plus() const { return sqrt(params.M*params.M-params.a*params.a)/(sqr(params.M+sqrt(params.M*params.M - params.a*params.a))+ sqr(params.a) ); }
Real KerrMetric::kappa_minus() const { return sqrt(params.M*params.M-params.a*params.a)/(sqr(params.M-sqrt(params.M*params.M - params.a*params.a))+ sqr(params.a) ) ; }
// conformal rescaling parameter getters
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