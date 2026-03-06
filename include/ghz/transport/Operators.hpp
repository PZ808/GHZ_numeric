//
// Created by Peter Zimmerman on 05.03.26.
//

#ifndef HS1DATA_DATA_OPERATORS_HPP
#define HS1DATA_DATA_OPERATORS_HPP

#pragma once
#include <vector>
#include <map>
#include <string>

#include "Corrector.hpp"

namespace ghz {

    using StateVec = ode::StateVec;
    using StateMat = std::vector<std::vector<StateVec>>;
    using RVector  = std::vector<Real>;
    using GHPSpectral = spectral::SpectralGHPVectorized;

    struct ModesMK { int m=0, kr=0, kz=0; };
    struct GHPType { int p=0, q=0; };

    struct XDerivs {
        GHPSpectral X;          // (p,q)
        GHPSpectral P_X;        // thorn  : (p+1,q+1)   (your convention: P = ∂r)
        GHPSpectral PpHr_X;      // thorn' : (p-1,q-1)   (Held thorn' on X)
        GHPSpectral PpHr_P_X;    // thorn'(P_X): (p, q)  since P_X is (p+1,q+1) then thorn' lowers -> (p,q)

        GHPSpectral EH_X;       // edthH(X)    : (p,   q-2)
        GHPSpectral EbH_X;      // edthBarH(X) : (p-2, q)

        GHPSpectral EbH_EH_X;   // ebh(eh(X))  : (p-2, q-2)
        GHPSpectral EH_EbH_X;   // eh(ebh(X))  : (p-2, q-2)

        GHPSpectral EH_P_X;     // eh(P_X)     : (p+1, q-1)
        GHPSpectral EbH_P_X;    // ebh(P_X)    : (p-1, q+1)

        XDerivs(size_t Nr, size_t Nz, ModesMK modes, GHPType type)
                : X       (Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p,   type.q)),
                  P_X     (Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p+1, type.q+1)),
                  PpHr_X   (Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p-1, type.q-1)),
                  PpHr_P_X (Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p,   type.q)),

                  EH_X    (Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p,   type.q-2)),
                  EbH_X   (Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p-2, type.q)),

                  EbH_EH_X(Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p-2, type.q-2)),
                  EH_EbH_X(Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p-2, type.q-2)),

                  EH_P_X  (Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p+1, type.q-1)),
                  EbH_P_X (Nr, Nz, {modes.m, modes.kr, modes.kz}, GHPScalar<Complex>(teuk::zeroC, type.p-1, type.q+1))
        {}
    };

    class GHPDifferX {
    public:

        GHPDifferX(const spectral::SpectralDiffer& diff,
                   const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords>& held,
                   const KinnersleyTetrad<OutgoingCoords>& tetrad)
                : diff_(diff), held_(held), tetrad_(tetrad),
                  Nz_(diff_.lgl_nodes().size()), Nr_(diff_.cl_nodes().size()) {}

        XDerivs diff_solution(const StateMat& state_2D,
                      const RVector& r_grid,
                      const RVector & z_grid,
                      ModesMK modes,
                      GHPType type,
                      Layout4 layout = {}) const;

    private:
        const spectral::SpectralDiffer& diff_;
        const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords>& held_;
        const KinnersleyTetrad<OutgoingCoords>& tetrad_;
        const size_t Nz_{};
        const size_t Nr_{};
    };

} // namespace ghz
#endif //HS1DATA_DATA_OPERATORS_HPP
