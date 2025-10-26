//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/KinnersleyTetrad.hpp"
#include "../include/SpinCoeffsNP.hpp"
#include "../include/MathMacros.hpp"
#include <complex>
#include <cmath>

using namespace math;

void KinnersleyTetradBL::build(double time, double r, double theta, double phi_azi) {
    (void)time; (void)phi_azi; // Kinnersley tetrad depends only on r, theta for Kerr

    std::complex<double> I(0.0,1.0);

    Real a = metric.a();
    Real M = metric.M();
    Real del = metric.Delta(r);
    Real sig = metric.Sigma(r, theta);

    // l.n = 1, m.mbar = -1
    l = { (r*r + a*a)/del, 1.0, 0.0, a/del };
    n = { (r*r + a*a)/(2.0*sig), -del/(2.0*sig), 0.0, a/(2.0*sig) };
    auto common = 1.0 / (std::sqrt(2.0)*(r + I*a*std::cos(theta)));
    m = ghz::CVector4{ I*a*std::sin(theta), 0.0, 1.0, I/std::sin(theta) } * common;
    mbar = m.conj();

    // Analytic spin coefficients
    Complex rho = -1.0/(r - I*a*std::cos(theta));
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc.set(SpinCoeffType::rho, rho);
    sc.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0);
    sc.set(SpinCoeffType::tau, -I*a*std::sin(theta)*rho*rhobar/std::sqrt(2.0));
    sc.set(SpinCoeffType::pi, I*a*std::sin(theta)*rho*rho/std::sqrt(2.0));
    sc.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0);
    sc.set(SpinCoeffType::beta, -rhobar*std::cos(theta)/(2.0*std::sin(theta)*std::sqrt(2.0)) );
    sc.set(SpinCoeffType::alpha, -std::conj(sc.get(SpinCoeffType::beta)) + sc.get(SpinCoeffType::pi));

    sc.set(SpinCoeffType::kappa, 0.0);
    sc.set(SpinCoeffType::sigma, 0.0);
    sc.set(SpinCoeffType::lambda, 0.0);
    sc.set(SpinCoeffType::nu, 0.0);
    sc.set(SpinCoeffType::epsilon, 0.0);

    //SpinCoefficientsGHP sc_ghp(sc);
    sc_ghp = SpinCoefficientsGHP(sc);
    sc_held = HeldCoefficients(sc_ghp, weyls);

    weyls.set(WeylScalarType::Psi0,0.0);
    weyls.set(WeylScalarType::Psi1,0.0);
    weyls.set(WeylScalarType::Psi2,M*cube(rho));
    weyls.set(WeylScalarType::Psi3,0.0);
    weyls.set(WeylScalarType::Psi4,0.0);
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
    m = ghz::CVector4{ I*a*std::sin(theta), 0.0, 1.0, I/std::sin(theta) } * common;
    mbar = m.conj();
}