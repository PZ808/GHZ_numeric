//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/KinnersleyTetrad.hpp"
#include "../include/SpinCoeffsNP.hpp"
#include "../include/MathMacros.hpp"
#include <complex>
#include <cmath>

using namespace math;
using namespace teuk::literals;

template <>
Complex KinnersleyTetrad<OutgoingCoords>::rho_at(const Real&  r, const Real& z) const {
    return -1.0_r/(r-I*z);
}
template <>
Complex KinnersleyTetrad<IngoingCoords>::rho_at(const Real&  r, const Real& z) const {
    return -1.0_r/(r-I*z);
}
template <>
Complex KinnersleyTetrad<BLCoords>::rho_at(const Real&  r, const Real& z) const {
    return -1.0_r/(r-I*z);
}

template <>
void KinnersleyTetrad<BLCoords>::build_tetrad_at(const BLCoords& Xbl) {

    using teuk::I;
    using namespace teuk::literals;

    // save cached coords for safety and consistency
    X_cached_ = Xbl;
    cache_valid_ = true;

    Real r = Xbl.x1;
    Real z = Xbl.x2; // z = \cos\theta

    Real a = metric.a();
    Real M = metric.M();
    Real del = metric.Delta(r);
    Real sig = metric.Sigma_z(r,z);
    Real s1 = math::Sqrt(1.0_r - z*z); // sin(th)
    Real s2 = 1.0_r - z*z; // sin^2(th)

    // ... compute tetrad members ...
    // l.n = 1, m.mbar = -1
    l = { (r*r + a*a)/del, 1.0_r, 0.0_r, a/del }; // contravariant vector
    n = { (r*r + a*a)/(2.0_r*sig), -del/(2.0_r*sig), 0.0_r, a/(2.0_r*sig) };
    auto common = 1.0_r / (sqrt(2.0_r)*(r + I*a*z));

    m = teuk::CVector4{I*a*s1, 0.0_r, -s1, I/s1} * common;
    mbar = m.conj();

    // Analytic spin coefficients
    Complex rho = -1.0_r/(r-I*a*z);
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc.set(SpinCoeffType::rho, rho);
    sc.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r);
    sc.set(SpinCoeffType::tau, -I*a*s1*rho*rhobar/math::Sqrt(2.0_r));
    sc.set(SpinCoeffType::pi, I*a*s1*rho*rho/math::Sqrt(2.0_r));
    sc.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0_r);
    sc.set(SpinCoeffType::beta, -rhobar*z/(two*s1*math::Sqrt(two)) );
    sc.set(SpinCoeffType::alpha, -std::conj(sc.get(SpinCoeffType::beta))
                                            + sc.get(SpinCoeffType::pi));

    sc.set(SpinCoeffType::kappa, teuk::zeroC);
    sc.set(SpinCoeffType::sigma, teuk::zeroC);
    sc.set(SpinCoeffType::lambda, teuk::zeroC);
    sc.set(SpinCoeffType::nu, teuk::zeroC);
    sc.set(SpinCoeffType::epsilon, teuk::zeroC);

    //SpinCoefficientsGHP sc_ghp(sc);
    sc_ghp = SpinCoefficientsGHP(sc);
    sc_held = HeldCoefficients(sc_ghp, weyls);

    weyls.set(WeylScalarType::Psi0,teuk::zeroC);
    weyls.set(WeylScalarType::Psi1,teuk::zeroC);
    weyls.set(WeylScalarType::Psi2,M*cube(rho));
    weyls.set(WeylScalarType::Psi3,teuk::zeroC);
    weyls.set(WeylScalarType::Psi4,teuk::zeroC);
}

template <>
void KinnersleyTetrad<OutgoingCoords>::build_tetrad_at(const OutgoingCoords & Xout) {

    using teuk::I;
    using namespace teuk::literals;

    // save cached coords for safety and consistency
    X_cached_ = Xout;
    cache_valid_ = true;


    Real a = metric.a();
    Real M = metric.M();

    Real r = Xout.x1;
    Real z = Xout.x2;

    Real del = metric.Delta(r);
    Real sig = sqr(a*z)+sqr(r);
    Real s2 = 1.0_r - sqr(z); // sin^2(th)
    Real s1 = sqrt(1.0_r-z*z); // sin(th)

    // contravariant NP null basis l.n = 1, m.mbar = -1
    l = { 0,1,0,0};
    n = { (r*r + a*a)/sig, -del/(2.0*sig), 0.0, a/sig };

    auto common = 1.0_r / (sqrt(2.0_r)*(r + I*a*z));
    m = teuk::CVector4{I * a * s1, 0.0_r, -s1 , I / s1 } * common;
    mbar = m.conj();

    Complex rho = -1.0_r/(r-I*a*z);
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc.set(SpinCoeffType::rho, rho);
    sc.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r);
    sc.set(SpinCoeffType::tau, -I*a*s1/(math::Sqrt(2.0_r)*sig));
    sc.set(SpinCoeffType::pi, I*a*s1*sqr(rho)/math::Sqrt(2.0_r));
    sc.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0_r);
    sc.set(SpinCoeffType::beta, -rhobar*z/(two*s1*math::Sqrt(two)) );
    sc.set(SpinCoeffType::alpha, sc.get(SpinCoeffType::beta) + sc.get(SpinCoeffType::pi));

    sc.set(SpinCoeffType::kappa, teuk::zeroC);
    sc.set(SpinCoeffType::sigma, teuk::zeroC);
    sc.set(SpinCoeffType::lambda, teuk::zeroC);
    sc.set(SpinCoeffType::nu, teuk::zeroC);
    sc.set(SpinCoeffType::epsilon, teuk::zeroC);

    //SpinCoefficientsGHP sc_ghp(sc);
    sc_ghp = SpinCoefficientsGHP(sc);
    sc_held = HeldCoefficients(sc_ghp, weyls);

    weyls.set(WeylScalarType::Psi0,0.0_r);
    weyls.set(WeylScalarType::Psi1,0.0_r);
    weyls.set(WeylScalarType::Psi2,M*cube(rho));
    weyls.set(WeylScalarType::Psi3,0.0_r);
    weyls.set(WeylScalarType::Psi4,0.0_r);
}

template <>
void KinnersleyTetrad<OutgoingCoordsCompact>::build_tetrad_at(const OutgoingCoordsCompact &Xout_C) {

    using teuk::I;
    using namespace teuk::literals;

    // save cached coords for safety and consistency
    X_cached_ = Xout_C;
    cache_valid_ = true;


    Real a = metric.a();
    Real M = metric.M();
    Real r = coord_helper.r_from_sigma(Xout_C.x1);  // r(\sigma) = \lambda/\sigma
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
    Complex delta_Ups = zero;

    // contravariant NP null basis l.n = 1, m.mbar = -1
    l = {0,1,0,0} ;
    n = { (r*r + a*a)/sig, -del/(2.0*sig)*dsigma_dr, 0.0, a/sig };

    // conformally rescale
    l = l/ sqr(Om_C);
    n = n * sqr(Om_C); //

    auto common = 1.0_r / (math::Sqrt(2.0_r)*(r + I*a*z));
    m = teuk::CVector4{I * a * s1, teuk::zeroC, -s1 , I / s1 } * common;
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

    sc.set(SpinCoeffType::kappa, 0.0_r);
    sc.set(SpinCoeffType::sigma, 0.0_r);
    sc.set(SpinCoeffType::lambda, 0.0_r);
    sc.set(SpinCoeffType::nu, 0.0_r);
    sc.set(SpinCoeffType::epsilon, 0.0_r);

    //SpinCoefficientsGHP sc_ghp(sc);
    sc_ghp = SpinCoefficientsGHP(sc);
    sc_held = HeldCoefficients(sc_ghp, weyls);

    weyls.set(WeylScalarType::Psi0,0.0_r);
    weyls.set(WeylScalarType::Psi1,0.0_r);
    weyls.set(WeylScalarType::Psi2,M*cube(rho)/cube(Om_C));
    weyls.set(WeylScalarType::Psi3,0.0_r);
    weyls.set(WeylScalarType::Psi4,0.0_r);
}

template <>
Tetrad::Scalars KinnersleyTetrad<BLCoords>::get_scalars_at(
        const BLCoords& Xbl) const {

    assert(cache_valid_ && "KinnersleyTetrad: build() must be called first");
    assert(X_cached_.has_value());
    assert(*X_cached_ == Xbl && "Mismatch: tetrad cached at different coordinates");

    using namespace teuk::literals;

    Real r = Xbl.x1;
    Real theta = Xbl.x2;

    Real a = metric.a();
    Real M = metric.M();
    Real del = metric.Delta(r);
    Real sig = metric.Sigma(r, theta);

    // --- Build the coefficients exactly as in build_tetrad() ---
    SpinCoefficients sc_local;
    // Analytic spin coefficients
    Complex rho = -1.0_r / (r - I * a * cos(theta));
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc_local.set(SpinCoeffType::rho, rho);
    sc_local.set(SpinCoeffType::mu, rhobar * sqr(rho) * del / 2.0_r);
    sc_local.set(SpinCoeffType::tau, -I * a * math::Sin(theta) * rho * rhobar / math::Sqrt(2.0_r));
    sc_local.set(SpinCoeffType::pi, I * a * math::Sin(theta) * rho * rho / math::Sqrt(2.0_r));
    sc_local.set(SpinCoeffType::kappa, teuk::zeroC);
    sc_local.set(SpinCoeffType::sigma, teuk::zeroC);
    sc_local.set(SpinCoeffType::lambda, teuk::zeroC);
    sc_local.set(SpinCoeffType::nu, teuk::zeroC);
    sc_local.set(SpinCoeffType::epsilon, teuk::zeroC);
    sc_local.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu) + rho * rhobar * (r - M) / 2.0_r);
    sc_local.set(SpinCoeffType::beta, -rhobar * math::Cos(theta) / (two * math::Sin(theta) * math::Sqrt(two)));
    sc_local.set(SpinCoeffType::alpha, -std::conj(sc.get(SpinCoeffType::beta)) + sc.get(SpinCoeffType::pi));

    //SpinCoefficientsGHP sc_ghp(sc);
    auto sc_ghp_local = SpinCoefficientsGHP(sc_local);
    WeylScalars W_local;

    W_local.set(WeylScalarType::Psi0, teuk::zeroC);
    W_local.set(WeylScalarType::Psi1, teuk::zeroC);
    W_local.set(WeylScalarType::Psi2, M * cube(rho));
    W_local.set(WeylScalarType::Psi3, teuk::zeroC);
    W_local.set(WeylScalarType::Psi4, teuk::zeroC);
    auto sc_held_local = HeldCoefficients(sc_ghp_local, W_local);

    Tetrad::Scalars scalars;
    scalars.ghp_scalars = sc_ghp_local;
    scalars.held_scalars = sc_held_local;

    return scalars;
}

template <>
Tetrad::Scalars KinnersleyTetrad<OutgoingCoords>::get_scalars_at(
        const OutgoingCoords& Xout) const
{

    assert(cache_valid_ && "KinnersleyTetrad: build() must be called first");
    assert(X_cached_.has_value());
    if (!(*X_cached_ == Xout)) {
        std::cerr << "Mismatch in coordinates:\n";
        std::cerr << "Cached: " << X_cached_->x0 << ", " << X_cached_->x1
                  << ", " << X_cached_->x2 << ", " << X_cached_->x3 << "\n";
        std::cerr << "Current: " << Xout.x0 << ", " << Xout.x1
                  << ", " << Xout.x2 << ", " << Xout.x3 << "\n";
        assert(false && "Mismatch: tetrad cached at different coordinates");
    }

 //   assert(*X_cached_ == Xout && "Mismatch: tetrad cached at different coordinates");

    using namespace teuk::literals;


    Tetrad::Scalars scalars;
    Real r = Xout.x1;
    Real z = Xout.x2;
    Real del = metric.Delta(r);
    using teuk::I;
    using namespace teuk::literals;
    Real a = metric.a();
    Real M = metric.M();
    Real sig = sqr(a*z)+sqr(r);
    Real s2 = 1.0_r - sqr(z); // sin^2(th)
    Real s1 = sqrt(1.0_r-z*z); // sin(th)

    // --- Build the coefficients exactly as in build_tetrad() ---
    SpinCoefficients sc_local;

    Complex rho = -1.0_r/(r-I*a*z);
    Complex rhobar = std::conj(rho);

    // Kinnersley tetrad BL values taken from Teukolsky Eq. (4.5)
    sc_local.set(SpinCoeffType::rho, rho);
    sc_local.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r);
    sc_local.set(SpinCoeffType::tau, -I*a*s1/(math::Sqrt(2.0_r)*sig));
    sc_local.set(SpinCoeffType::pi, I*a*s1*sqr(rho)/math::Sqrt(2.0_r));
    sc_local.set(SpinCoeffType::gamma, sc.get(SpinCoeffType::mu)+ rho*rhobar*(r-M)/2.0_r);
    sc_local.set(SpinCoeffType::beta, -rhobar*z/(two*s1*math::Sqrt(two)) );
    sc_local.set(SpinCoeffType::alpha, sc.get(SpinCoeffType::beta) + sc.get(SpinCoeffType::pi));

    sc_local.set(SpinCoeffType::kappa, teuk::zeroC);
    sc_local.set(SpinCoeffType::sigma, teuk::zeroC);
    sc_local.set(SpinCoeffType::lambda, teuk::zeroC);
    sc_local.set(SpinCoeffType::nu, teuk::zeroC);
    sc_local.set(SpinCoeffType::epsilon, teuk::zeroC);

    WeylScalars weyls_local;
    weyls_local.set(WeylScalarType::Psi0,0.0_r);
    weyls_local.set(WeylScalarType::Psi1,0.0_r);
    weyls_local.set(WeylScalarType::Psi2,M*cube(rho));
    weyls_local.set(WeylScalarType::Psi3,0.0_r);
    weyls_local.set(WeylScalarType::Psi4,0.0_r);

    auto sc_ghp_local = SpinCoefficientsGHP(sc_local);
    HeldCoefficients sc_held_local = HeldCoefficients(sc_ghp_local, weyls_local);
    scalars.ghp_scalars = sc_ghp_local;
    scalars.held_scalars = sc_held_local;

    return scalars;
}

template <>
Tetrad::Scalars KinnersleyTetrad<OutgoingCoordsCompact>::get_scalars_at(
        const OutgoingCoordsCompact& XoutC) const
{
    using namespace teuk::literals;
    Tetrad::Scalars scalars;

    Real a  = metric.a();
    Real M  = metric.M();
    Real r  = coord_helper.r_from_sigma(XoutC.x1);
    Real z  = XoutC.x2;

    Real del = metric.Delta(r);
    Real sig = sqr(a*z)+sqr(r);
    Real Om_C = XoutC.x1/metric.lambda_C();
    Real s1 = sqrt(1.0_r-z*z);

    Real rho0_C = (metric.r_plus()/metric.lambda_C());
    Real dsigma_dr = -metric.lambda_C()*rho0_C/(r*r);
    Real dOm_dr = 1.0_r/metric.lambda_C()*dsigma_dr;

    Real Delta_Ups = - del/(2.0*sig)*sqr(Om_C);
    Complex rho = -1.0_r/(r - I*a*z);
    Complex rhobar = std::conj(rho);

    // --- Build the coefficients exactly as in build_tetrad() ---
    SpinCoefficients sc_local;

    sc_local.set(SpinCoeffType::rho, rho/sqr(Om_C));
    sc_local.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r + dOm_dr/Om_C);
    sc_local.set(SpinCoeffType::tau, -I*a*s1/(sqrt(2)*sig)/Om_C);
    sc_local.set(SpinCoeffType::pi, I*a*s1*sqr(rho)/sqrt(2.0_r)/Om_C);

    sc_local.set(SpinCoeffType::gamma,
                 sc_local.get(SpinCoeffType::mu) + rho*rhobar*(r-M)/2.0_r - Delta_Ups);

    sc_local.set(SpinCoeffType::beta, -(rhobar*z/(2.0_r*s1*sqrt(2.0_r)))/Om_C);
    sc_local.set(SpinCoeffType::alpha,
                 (sc_local.get(SpinCoeffType::beta)
                  + sc_local.get(SpinCoeffType::pi)) / Om_C);

    sc_local.set(SpinCoeffType::kappa, 0.0_r);
    sc_local.set(SpinCoeffType::sigma, 0.0_r);
    sc_local.set(SpinCoeffType::lambda, 0.0_r);
    sc_local.set(SpinCoeffType::nu, 0.0_r);
    sc_local.set(SpinCoeffType::epsilon, 0.0_r);

    WeylScalars W;
    W.set(WeylScalarType::Psi0, 0.0_r);
    W.set(WeylScalarType::Psi1, 0.0_r);
    W.set(WeylScalarType::Psi2, M*cube(-1.0_r/(r))*cube(1.0/Om_C)); // same as your code
    W.set(WeylScalarType::Psi3, 0.0_r);
    W.set(WeylScalarType::Psi4, 0.0_r);
    auto sc_ghp_local = SpinCoefficientsGHP(sc_local);
    HeldCoefficients sc_held_local = HeldCoefficients(sc_ghp_local, W);
    scalars.ghp_scalars = sc_ghp_local;
    scalars.held_scalars = sc_held;
    return scalars;
}

template <>
SpinCoefficientsGHP KinnersleyTetrad<OutgoingCoordsCompact>::get_ghp_scalars_at(
        const OutgoingCoordsCompact& X) const
{
    using namespace teuk::literals;

    Real a  = metric.a();
    Real M  = metric.M();
    Real r  = coord_helper.r_from_sigma(X.x1);
    Real z  = X.x2;

    Real del = metric.Delta(r);
    Real sig = sqr(a*z)+sqr(r);
    Real Om_C = X.x1/metric.lambda_C();
    Real s1 = sqrt(1.0-z*z);

    Real rho0_C = (metric.r_plus()/metric.lambda_C());
    Real dsigma_dr = -metric.lambda_C()*rho0_C/(r*r);
    Real dOm_dr = 1.0/metric.lambda_C()*dsigma_dr;

    Real Delta_Ups = - del/(2.0*sig)*sqr(Om_C);
    Complex rho = -1.0_r/(r - I*a*z);
    Complex rhobar = std::conj(rho);

    // --- Build the coefficients exactly as in build_tetrad() ---
    SpinCoefficients sc_local;

    sc_local.set(SpinCoeffType::rho, rho/sqr(Om_C));
    sc_local.set(SpinCoeffType::mu, rhobar*sqr(rho)*del/2.0_r + dOm_dr/Om_C);
    sc_local.set(SpinCoeffType::tau, -I*a*s1/(sqrt(2)*sig)/Om_C);
    sc_local.set(SpinCoeffType::pi, I*a*s1*sqr(rho)/sqrt(2.0_r)/Om_C);

    sc_local.set(SpinCoeffType::gamma,
                 sc_local.get(SpinCoeffType::mu) + rho*rhobar*(r-M)/2.0_r - Delta_Ups);

    sc_local.set(SpinCoeffType::beta, -(rhobar*z/(2.0*s1*sqrt(2.0)))/Om_C);
    sc_local.set(SpinCoeffType::alpha,
                 (sc_local.get(SpinCoeffType::beta)
                  + sc_local.get(SpinCoeffType::pi)) / Om_C);

    sc_local.set(SpinCoeffType::kappa, 0.0_r);
    sc_local.set(SpinCoeffType::sigma, 0.0_r);
    sc_local.set(SpinCoeffType::lambda, 0.0_r);
    sc_local.set(SpinCoeffType::nu, 0.0_r);
    sc_local.set(SpinCoeffType::epsilon, 0.0_r);

    return SpinCoefficientsGHP(sc_local); // generate GHP coefficients from NP ones
}

template <>
WeylScalars
KinnersleyTetrad<OutgoingCoordsCompact>::get_weyl_scalars_at(
        const OutgoingCoordsCompact& X) const
{
    using namespace teuk::literals;

    Real r = coord_helper.r_from_sigma(X.x1);
    Real Om_C = X.x1 / metric.lambda_C();
    Real M = metric.M();

    WeylScalars W;
    W.set(WeylScalarType::Psi0, 0.0_r);
    W.set(WeylScalarType::Psi1, 0.0_r);
    W.set(WeylScalarType::Psi2, M*cube(-1.0_r/(r))*cube(1.0/Om_C)); // same as your code
    W.set(WeylScalarType::Psi3, 0.0_r);
    W.set(WeylScalarType::Psi4, 0.0_r);

    return W;
}

/**
 * @brief Compute the full Held coefficients at coordinates X.
 *
 * This templated implementation works for *any* coordinate chart for which the
 * tetrad provides:
 *    - get_spin_coeffs_at(X)
 *    - get_weyl_scalars_at(X)
 *
 * This avoids redundant explicit template specializations for each coordinate
 * system (BLCoords, OutgoingCoords, OutgoingCoordsCompact, etc.).
 *
 * @tparam CoordT  The coordinate type
 * @param  X       Coordinates where the coefficients are evaluated
 * @return HeldCoefficients  Combined spin coefficients + Weyl scalars
*/
template <typename CoordT>
HeldCoefficients KinnersleyTetrad<CoordT>::get_held_scalars_at(const CoordT &X) const {

    return HeldCoefficients(get_spin_coeffs_at(X),
                            get_weyl_scalars_at(X));
}

// for any CoordT you actually use:
//template class KinnersleyTetrad<OutgoingCoordsCompact>;
