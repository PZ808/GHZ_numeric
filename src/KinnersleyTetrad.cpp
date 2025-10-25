//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/KinnersleyTetrad.hpp"
#include "../include/SpinCoeffsNP.hpp"
#include <complex>
#include <cmath>

void KinnersleyTetradBL::build(double time, double r, double theta, double phi_azi) {
    (void)time; (void)phi_azi; // Kinnersley tetrad depends only on r, theta for Kerr

    double a = metric.a();
    double del = metric.Delta(r);
    double sig = metric.Sigma(r, theta);

    l = { (r*r + a*a)/del, 1.0, 0.0, a/del };
    n = { (r*r + a*a)/(2.0*sig), -del/(2.0*sig), 0.0, a/(2.0*sig) };

    std::complex<double> I(0.0,1.0);
    auto common = 1.0 / (std::sqrt(2.0)*(r + I*a*std::cos(theta)));
    m = CVector4{ I*a*std::sin(theta), 0.0, 1.0, I/std::sin(theta) } * common;
    mbar = m.conj();

    // Analytic spin coefficients
    Complex rho = -1.0/(r - I*a*std::cos(theta));
    Complex rhobar = std::conj(rho);

    sc.set(SpinCoeffType::rho, rho);
    sc.set(SpinCoeffType::mu, rhobar*rhobar*rho*del/(2.0*sig));
    sc.set(SpinCoeffType::gamma, rho*rhobar*metric.M()/2.0);
    sc.set(SpinCoeffType::tau, -I*a*std::sin(theta)*rho*rho/std::sqrt(2.0));
    sc.set(SpinCoeffType::pi, I*a*std::sin(theta)*rho*rhobar/std::sqrt(2.0));
    sc.set(SpinCoeffType::beta, -rho*std::tan(theta/2.0)/(2*std::sqrt(2.0)));
    sc.set(SpinCoeffType::alpha, std::conj(sc.get(SpinCoeffType::beta)) - sc.get(SpinCoeffType::pi));

    sc.set(SpinCoeffType::kappa, 0.0);
    sc.set(SpinCoeffType::sigma, 0.0);
    sc.set(SpinCoeffType::lambda, 0.0);
    sc.set(SpinCoeffType::nu, 0.0);
    sc.set(SpinCoeffType::epsilon, 0.0);
}


void KinnersleyTetradOutgoing::build(double time, double r, double theta, double phi_azi) {
    (void)time; (void)phi_azi; // Kinnersley tetrad depends only on r, theta for Kerr

    double a = metric.a();
    double del = metric.Delta(r);
    double sig = metric.Sigma(r, theta);

    l = { (r*r + a*a)/del, 1.0, 0.0, a/del };
    n = { (r*r + a*a)/(2.0*sig), -del/(2.0*sig), 0.0, a/(2.0*sig) };

    std::complex<double> I(0.0,1.0);
    auto common = 1.0 / (std::sqrt(2.0)*(r + I*a*std::cos(theta)));
    m = CVector4{ I*a*std::sin(theta), 0.0, 1.0, I/std::sin(theta) } * common;
    mbar = m.conj();
}