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
    using namespace teuk::literals;
    using teuk::I;

    Real a = metric.a();
    Real M = metric.M();
    Real del = metric.Delta(r);
    Real sig = metric.Sigma(r, theta);

    // l.n = 1, m.mbar = -1
    l = { (r*r + a*a)/del, 1.0, 0.0, a/del }; // contravariant vector
    n = { (r*r + a*a)/(2.0*sig), -del/(2.0*sig), 0.0, a/(2.0*sig) };
    auto common =1.0_r / (math::sqrt(Real(2.0))*(r + I*a*math::cos(theta)));
    m = ghz::CVector4{ I*a*math::sin(theta), 0.0, 1.0, I/math::sin(theta) } * common;
    mbar = m.conj();

    // Analytic spin coefficients
    Complex rho = -Real(1.0)/(r - I*a*math::cos(theta));
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc.set(SpinCoeffType::rho, rho);
    sc.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r);
    sc.set(SpinCoeffType::tau, -I*a*sin(theta)*rho*rhobar/sqrt(2.0_r));
    sc.set(SpinCoeffType::pi, I*a*sin(theta)*rho*rho/sqrt(2.0_r));
    sc.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0_r);
    sc.set(SpinCoeffType::beta, -rhobar*cos(theta)/(2.0*sin(theta)*sqrt(2.0_r)) );
    sc.set(SpinCoeffType::alpha, -std::conj(sc.get(SpinCoeffType::beta)) + sc.get(SpinCoeffType::pi));

    sc.set(SpinCoeffType::kappa, 0.0_r);
    sc.set(SpinCoeffType::sigma, 0.0_r);
    sc.set(SpinCoeffType::lambda, 0.0_r);
    sc.set(SpinCoeffType::nu, 0.0_r);
    sc.set(SpinCoeffType::epsilon, 0.0_r);

    //SpinCoefficientsGHP sc_ghp(sc);
    sc_ghp = SpinCoefficientsGHP(sc);
    sc_held = HeldCoefficients(sc_ghp, weyls);

    weyls.set(WeylScalarType::Psi0,0.0);
    weyls.set(WeylScalarType::Psi1,0.0);
    weyls.set(WeylScalarType::Psi2,M*cube(rho));
    weyls.set(WeylScalarType::Psi3,0.0);
    weyls.set(WeylScalarType::Psi4,0.0);
}


void KinnersleyTetradBL::build_tetrad(const BLCoords& Xbl) {

    using namespace teuk::literals;

    Real r = Xbl.x1;
    Real theta = Xbl.x2;

    Real a = metric.a();
    Real M = metric.M();
    Real del = metric.Delta(r);
    Real sig = metric.Sigma(r, theta);

    // l.n = 1, m.mbar = -1
    l = { (r*r + a*a)/del, 1.0_r, 0.0_r, a/del }; // contravariant vector
    n = { (r*r + a*a)/(2.0_r*sig), -del/(2.0_r*sig), 0.0_r, a/(2.0_r*sig) };
    auto common = 1.0_r / (sqrt(2.0_r)*(r + I*a*cos(theta)));
    m = ghz::CVector4{ I*a*std::sin(theta), 0.0_r, 1.0_r, I/std::sin(theta) } * common;
    mbar = m.conj();

    // Analytic spin coefficients
    Complex rho = -1.0_r/(r - I*a*cos(theta));
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc.set(SpinCoeffType::rho, rho);
    sc.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r);
    sc.set(SpinCoeffType::tau, -I*a*std::sin(theta)*rho*rhobar/sqrt(2.0_r));
    sc.set(SpinCoeffType::pi, I*a*std::sin(theta)*rho*rho/sqrt(2.0_r));
    sc.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0_r);
    sc.set(SpinCoeffType::beta, -rhobar*cos(theta)/(2.0*sin(theta)*sqrt(2.0)) );
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

//void KinnersleyTetradBL::build_tetrad_compact(teuk::Real u, teuk::Real sigma, teuk::Real z, teuk::Real ph) {
    // r direction compactified
    // to be implemented
//}

    void KinnersleyTetradOutgoing::build(Real time, Real r, Real z, Real phi_azi) {
    (void)time; (void)phi_azi; // builds Kinnersley tetrad and derived quantities in  outgoing coordinates (u,r,z,phi)

    using teuk::I;
    using namespace teuk::literals;
    Real a = metric.a();
    Real M = metric.M();
    Real del = metric.Delta(r);
    Real sig = sqr(a*z)+sqr(r);
    Real s2 = 1.0_r - sqr(z); // sin^2(th)
    Real s1 = sqrt(1.0_r-z*z); // sin(th)


    // contravariant NP null basis l.n = 1, m.mbar = -1
    l = { 0,1,0,0};
    n = { (r*r + a*a)/sig, -del/(2.0*sig), 0.0, a/sig };


    auto common = 1.0_r / (sqrt(2.0_r)*(r + I*a*z));

    m = ghz::CVector4{ I*a*s1, 0.0, -s1 , I/s1 } * common;
    mbar = m.conj();


    Complex rho = -1.0_r/(r-I*a*z);
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc.set(SpinCoeffType::rho, rho);
    sc.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r);
    sc.set(SpinCoeffType::tau, -I*a*s1/(sqrt(2.0_r)*sig));
    sc.set(SpinCoeffType::pi, I*a*s1*sqr(rho)/sqrt(2.0_r));
    sc.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0_r);
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

void KinnersleyTetradOutgoing::build_tetrad_compact(const OutgoingCoordsCompact &Xout_C) {
    // sigma = \lambda rho0/r
    using teuk::I;
    using namespace teuk::literals;

    Real a = metric.a();
    Real M = metric.M();
    Real r = coords.r_from_sigma(Xout_C.x1);  // r(\sigma) = \lambda/\sigma
    Real z = Xout_C.x2;
    Real del = metric.Delta(r);
    Real sig = sqr(a*z)+sqr(r);
    Real s2 = 1.0 - sqr(z); // sin^2(th)
    Real Om_C = Xout_C.x1/metric.lambda_C(); // choice of conformal factor \Omega = \sigma/\lambda = 1/r
    Real s1 = sqrt(1.0-z*z); // sin(th)
    Real rho0_C = (metric.r_plus()/metric.lambda_C());;
    Real dsigma_dr = -metric.lambda_C()*rho0_C / (r*r);
    Real dOm_dr = 1.0/metric.lambda_C()*dsigma_dr;
    Real Ups = - Om_C; // Upsilon = Om_C^{-1} d\Omega_C/dr
    Real Delta_Ups = - del/(2.0*sig)*sqr(Om_C);
    Complex delta_Ups = 0;



    // contravariant NP null basis l.n = 1, m.mbar = -1
    l = {0,1,0,0} ;
    n = { (r*r + a*a)/sig, -del/(2.0*sig)*dsigma_dr, 0.0, a/sig };

    // conformally rescale
    l = l/ sqr(Om_C);
    n = n * sqr(Om_C); //

    auto common = 1.0_r / (std::sqrt(2.0_r)*(r + I*a*z));
    m = ghz::CVector4{ I*a*s1, 0.0, -s1 , I/s1 } * common;
    mbar = m.conj();


    Complex rho = -1.0_r/(r-I*a*z);
    Complex rhobar = std::conj(rho);

    // transform values take from Steward's Advanced GR book (3.11.2) pp 156
    sc.set(SpinCoeffType::rho, rho/sqr(Om_C)) ;
    sc.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r+ dOm_dr/Om_C);
    sc.set(SpinCoeffType::tau, -I*a*s1/(sqrt(2)*sig)/Om_C);
    sc.set(SpinCoeffType::pi, I*a*s1*sqr(rho)/sqrt(2.0_r)/Om_C + delta_Ups);
    sc.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0_r - Delta_Ups);
    sc.set(SpinCoeffType::beta, -(rhobar*z/(2.0*s1*sqrt(2.0))) / Om_C);
    sc.set(SpinCoeffType::alpha, (sc.get(SpinCoeffType::beta) + sc.get(SpinCoeffType::pi) - conj(delta_Ups)) / Om_C);

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
    weyls.set(WeylScalarType::Psi2,M*cube(rho)/cube(Om_C));
    weyls.set(WeylScalarType::Psi3,0.0);
    weyls.set(WeylScalarType::Psi4,0.0);


}