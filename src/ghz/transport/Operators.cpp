//
// Created by Peter Zimmerman on 05.03.26.
//

#include "../include/ghz/transport/Operators.hpp"

#include "../include/ghz/spectral/KinnersleySpectralHeldOperators.hpp"   // KinnersleyHeldOperators
#include "../include/ghz/ghp/GHPScalars.hpp"

#include <cassert>
#include <omp.h>

    namespace ghz {

        XDerivs GHPDifferX::diff_solution(
                const StateMat& state_2D,
                const RVector& r_grid,
                const RVector& z_grid,
                ModesMK modes,
                GHPType type,
                Layout4 L) const
        {
            const size_t Nr = r_grid.size();
            const size_t Nz = z_grid.size();

            // Sanity: state_2D is indexed [iz][ir]
            assert(state_2D.size() == Nz);
            if (Nz > 0) assert(state_2D[0].size() == Nr);

            // Prefer a single source of truth for z nodes: if you want strict consistency with diff_,
            // assert equality (or just ignore z_grid and use diff_.lgl_nodes()).
            // Here: accept z_grid but guard mismatch.
            assert(diff_.lgl_nodes().size() == Nz);

            const int p  = type.p;
            const int q  = type.q;
            const int m  = modes.m;
            const int kr = modes.kr;
            const int kz = modes.kz;

            // --- Allocate outputs (strongly typed container) ---
            XDerivs D;

            // Base field X has type (p,q)
            D.X = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p, q));

            // thorn X (thorn = ∂r) has type (p+1,q+1)
            D.P_X = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p+1, q+1));

            // thorn' Held (thornPH = ∂u-Delta/(2Sigma)∂r) has type (p-1,q-1)
            //D.PpH_X = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-1, q-1));
            D.PpH_X = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-1, q-1));

            // Held edth outputs:
            // edthH : (p,q) -> (p, q-2)
            // edthBarH : (p,q) -> (p-2, q)
            D.EH_X  = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p,   q-2));
            D.EbH_X = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-2, q));

            // Apply edth/edthBar to thorn_X (which is (p+1,q+1)):
            // edthH(thorn_X):    (p+1,q+1) -> (p+1, q-1)
            // edthBarH(thorn_X): (p+1,q+1) -> (p-1, q+1)
            D.EH_P_X  = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p+1, q-1));
            D.EbH_P_X = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-1, q+1));

            // Compositions on X:
            // edthH(edthBarH(X)): (p-2,q) -> (p-2, q-2)
            // edthBarH(edthH(X)): (p,q-2) -> (p-2, q-2)
            D.EH_EbH_X = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-2, q-2));
            D.EbH_EH_X = GHPSpectral(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-2, q-2));

            // --- Fill X and thorn_X directly from ODE state ---
#pragma omp parallel for collapse(2) default(none) shared(D, state_2D, L, Nr, Nz, p, q)
            for (size_t iz = 0; iz < Nz; ++iz) {
                for (size_t ir = 0; ir < Nr; ++ir) {
                    const auto& y = state_2D[iz][ir];

                    const Complex Xval (y[L.re],  y[L.im]);
                    const Complex dXdr (y[L.reR], y[L.imR]);   // stored radial derivative

                    D.X.set_index(ir, iz, GHPScalar<Complex>(Xval, p, q));
                    D.P_X.set_index(ir, iz, GHPScalar<Complex>(dXdr, p+1, q+1));
                }
            }

            // --- Apply Held operators on each r-slice ---
            // Construct operators once (thread-safe if diff_/held_/tetrad_ are read-only)
            KinnersleyHeldOperators<OutgoingCoords> held_ops(diff_, held_, tetrad_);

#pragma omp parallel for default(none) shared(p, q, D, held_ops, Nr)
            for (size_t ir = 0; ir < Nr; ++ir) {

                auto slice_X  = D.X.slice_R(ir);
                auto slice_P_X= D.P_X.slice_R(ir);

                auto slice_PpH_X= D.PpH_X.slice_R(ir);

                auto slice_EH_X= D.EH_X.slice_R(ir);
                auto slice_EbH_X= D.EbH_X.slice_R(ir);

                auto slice_Eh_P_X= D.EH_P_X.slice_R(ir);
                auto slice_EbH_P_X= D.EbH_P_X.slice_R(ir);

                auto slice_EH_EbH_X= D.EH_EbH_X.slice_R(ir);
                auto slice_EbH_EH_X= D.EbH_EH_X.slice_R(ir);

                // thornPH on X (frequency-domain time derivative)
                //held_ops.thornPH_inplace_RSliceV(slice_X, slice_PpH_X);
                held_ops.thornPHr_inplace_RSliceV(tetrad_, slice_X,
                                                  slice_P_X, slice_PpH_X);

                // edth/edthBar on X
                held_ops.edthH_inplace_RSliceV(slice_X, slice_EH_X);
                held_ops.edthBarH_inplace_RSliceV(slice_X, slice_EbH_X);

                // compositions
                held_ops.edthH_inplace_RSliceV(slice_EbH_X, slice_EH_EbH_X);
                held_ops.edthBarH_inplace_RSliceV(slice_Eh_P_X, slice_EbH_EH_X);

                // edth/edthBar on thorn_X
                held_ops.edthH_inplace_RSliceV(slice_P_X, slice_Eh_P_X);
                held_ops.edthBarH_inplace_RSliceV(slice_P_X, slice_EbH_P_X);

#ifndef NDEBUG
                // Quick weight sanity checks on this slice (assumes weights uniform across slice)
                assert(slice_X[0].p() == p && slice_X[0].q() == q);
                assert(slice_EH_X[0].p() == p && slice_EH_X[0].q() == q-2);
                assert(slice_EbH_X[0].p() == p-2 && slice_EbH_X[0].q() == q);
                assert(slice_EH_EbH_X[0].p() == p-2 && slice_EH_EbH_X[0].q() == q-2);
                assert(slice_EbH_EH_X[0].p() == p-2 && slice_EbH_EH_X[0].q() == q-2);
#endif
            }

            return D;
        }

    } // namespace ghz

}