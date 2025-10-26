//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/KinnersleyTetrad.hpp"
#include "../include/SpinCoeffsNP.hpp"
#include "../include/MathMacros.hpp"
#include <complex>
#include <cmath>

using namespace math;

void KinnersleyTetradBL::build(Real time, Real r, Real theta, Real phi_azi) {
    (void)time; (void)phi_azi; // builds Kinnersley tetrad and derived quantities in BL coordinates

    using teuk::I;

    Real a = metric.a();
    Real M = metric.M();
    Real del = metric.Delta(r);
    Real sig = metric.Sigma(r, theta);

    // l.n = 1, m.mbar = -1
    l = { (r*r + a*a)/del, 1.0, 0.0, a/del }; // contravariant vector
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


void KinnersleyTetradOutgoing::build(Real time, Real r, Real z, Real phi_azi) {
    (void)time; (void)phi_azi; // builds Kinnersley tetrad and derived quantities in  outgoing coordinates (u,r,z,phi)

    using teuk::I;
    Real a = metric.a();
    Real M = metric.M();
    Real del = metric.Delta(r);
    Real sig = sqr(a*z)+sqr(r);
    Real s2 = 1.0 - sqr(z); // sin^2(th)
    Real s1 = sqrt(1.0-z*z); // sin(th)

    // contravariant NP null basis l.n = 1, m.mbar = -1
    l = { 1,0,0,0};
    n = { (r*r + a*a)/sig, -del/(2.0*sig), 0.0, a/sig };

    auto common = 1.0 / (std::sqrt(2.0)*(r + I*a*z));
    m = ghz::CVector4{ I*a*s1, 0.0, -s1 , I/s1 } * common;
    mbar = m.conj();


    Complex rho = -1.0/(r-I*a*z);
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc.set(SpinCoeffType::rho, rho);
    sc.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0);
    sc.set(SpinCoeffType::tau, -I*a*s1/(sqrt(2)*sig));
    sc.set(SpinCoeffType::pi, I*a*s1*sqr(rho)/sqrt(2));
    sc.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0);
    sc.set(SpinCoeffType::beta, -rhobar*z/(2.0*s1*sqrt(2.0)) );
    sc.set(SpinCoeffType::alpha, sc.get(SpinCoeffType::beta) + sc.get(SpinCoeffType::pi));

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