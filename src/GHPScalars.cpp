//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/GHPScalars.hpp"

SpinCoefficientsGHP::SpinCoefficientsGHP(const SpinCoefficients &sc_np) {
    // initialize weights according to GHP convention (p,q)
    // (using Held’s sign conventions)
    using SCT = SpinCoeffType;

    kappa  = GHPScalar(sc_np.get(SCT::kappa), 3, 1);
    kappap = GHPScalar(-sc_np.get(SCT::nu), -3, -1);
    sigma  = GHPScalar(sc_np.get(SCT::sigma), 3, -1);
    sigmap = GHPScalar(-sc_np.get(SCT::lambda), -3, 1);
    rho    = GHPScalar(sc_np.get(SCT::rho), 1, 1);
    rhop   = GHPScalar(-sc_np.get(SCT::mu), -1, -1);
    tau    = GHPScalar(sc_np.get(SCT::tau), 1, -1);
    taup   = GHPScalar(sc_np.get(SCT::pi), -1, 1);
}
