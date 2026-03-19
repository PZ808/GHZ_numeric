//
// Created by Peter Zimmerman on 25.10.25.
//

#include "ghz/ghp/GHPScalars.hpp"


ghp::SpinCoefficientsGHP::SpinCoefficientsGHP(const SpinCoefficients &sc_np) {
    // initialize weights according to GHP convention (p,q)
    // (using Held’s sign conventions)
    using SCT = SpinCoeffType;

    // boost/spin covariant definite weight scalars
    kappa  = GHPScalar(sc_np.get(SCT::kappa), 3, 1);
    kappa_bar  = GHPScalar(std::conj(sc_np.get(SCT::kappa)), 1, 3);
    kappap = GHPScalar(-sc_np.get(SCT::nu), -3, -1);
    kappap_bar = GHPScalar(-std::conj(sc_np.get(SCT::nu)), -1, -3);
    sigma  = GHPScalar(sc_np.get(SCT::sigma), 3, -1);
    sigma_bar  = GHPScalar(std::conj(sc_np.get(SCT::sigma)), -1, 3);
    sigmap = GHPScalar(-sc_np.get(SCT::lambda), -3, 1);
    sigmap_bar = GHPScalar(std::conj(-sc_np.get(SCT::lambda)), 1, -3);
    rho    = GHPScalar(sc_np.get(SCT::rho), 1, 1);
    rho_bar    = GHPScalar(std::conj(sc_np.get(SCT::rho)), 1, 1);
    rhop   = GHPScalar(-sc_np.get(SCT::mu), -1, -1);
    rhop_bar   = GHPScalar(std::conj(-sc_np.get(SCT::mu)), -1, -1);
    tau    = GHPScalar(sc_np.get(SCT::tau), 1, -1);
    tau_bar    = GHPScalar(std::conj(sc_np.get(SCT::tau)), -1, 1);
    taup   = GHPScalar(-sc_np.get(SCT::pi), -1, 1);
    taup_bar   = GHPScalar(std::conj(-sc_np.get(SCT::pi)), 1, -1);

    // indefinite weight scalars
    beta = sc_np.get(SCT::beta);
    betap = -sc_np.get(SCT::alpha);
    epsilon = sc_np.get(SCT::epsilon);
    epsilonp = -sc_np.get(SCT::gamma);
}

void ghp::GHPCoefficients::set(GHPCoefficientType type, GHPScalar<Complex> value) {
    coeffs[type] = value;
}

ghp::GHPScalar<Complex> ghp::GHPCoefficients::get(GHPCoefficientType type) const {
    auto it = coeffs.find(type);
    if (it != coeffs.end()) return it->second;
    return GHPScalar(teuk::zeroC,0,0);
}