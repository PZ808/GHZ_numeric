//
// Created by Peter Zimmerman on 05.03.26.
//
// SourceBuilders.cpp
//
// Collapsed + cleaned implementation of:
//   source_xmmb, source_xnm, source_xnn
//
// Assumptions:
//  - XDerivs has members:
//      X, P_X, PpH_X, PpH_P_X, EH_X, EbH_X, EbH_EH_X, EH_EbH_X, EH_P_X, EbH_P_X
//  - GHPScalar<Complex> supports +,*, conj(), .value(), etc.
//  - ghp::GHPBackgroundFieldsVectorized provides rho(ir,iz) as a GHPScalar<Complex> (or similar).
//
#include "../include/ghz/transport/SourceBuilders.hpp"
#include "ghz/core/GhzTypes.hpp"


namespace ghz::source {

    using namespace teuk::literals;
    Complex I = teuk::I;

    using GHPSpectral = spectral::SpectralGHPVectorized;

// -----------------------------------------------------------------------------
// x_mmbar source: just T_ll if present, else zeros
// -----------------------------------------------------------------------------

    GHPSpectral source_xmmb(
            size_t Nr, size_t Nz, ModesMK modes, GHPType out_type,
            const GHPSpectral *Tll) {
        if (Tll) return *Tll;

        return GHPSpectral(
                Nr, Nz,
                {modes.m, modes.kr, modes.kz},
                GHPScalar<Complex>(teuk::zeroC, out_type.p, out_type.q));
    }

// -----------------------------------------------------------------------------
// x_nm source: N[x_mmbar] + T_lm
// -----------------------------------------------------------------------------

    GHPSpectral source_xnm(
            const XDerivs &xmmbar,
            const std::vector<Real> &r_grid,
            const std::vector<Real> &z_grid,
            const KerrMetricOutgoing &metric,
            const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords> &held,
            const ghp::GHPBackgroundFieldsVectorized &ghp,
            GHPType out_type,
            const GHPSpectral *Tlm) {
        const size_t Nr = xmmbar.X.Nr();
        const size_t Nz = xmmbar.X.Nz();

        // Mode labels: assume all deriv fields share same modes
        const auto modes = xmmbar.X.modes(); // adjust accessor if needed

        GHPSpectral NX_source(
                Nr, Nz, modes,
                GHPScalar<Complex>(teuk::zeroC, out_type.p, out_type.q));

        const Real a = metric.a();

#pragma omp parallel for collapse(2) default(none) \
    shared(NX_source, xmmbar, held, ghp, r_grid, z_grid, half, three, four, out_type, Tlm, a, Nr, Nz)
        for (size_t iz = 0; iz < Nz; ++iz) {
            const auto tau0 = held.tauH(iz); // τ^°(z), Held-circle

            for (size_t ir = 0; ir < Nr; ++ir) {
                // Background rho from your precomputed GHP background fields
                const auto rho_ghp = ghp.rho(ir, iz);
                const auto rhob_ghp = rho_ghp.conj();

                // Geometric rho(r,z) (needed by the specific N[x] formula you used)
                // rho = -1/(r - i a z), but here we need r,z. Prefer using your stored grids.
                // We assume xmmbar.X stores/knows its grids; otherwise replace with your grid accessors.
                const Real r = r_grid[ir];
                const Real z = z_grid[iz];

                const Complex rho = -Real(1) / (Complex(r, -a * z)); // r - i a z
                const Complex rhob = std::conj(rho);

                // Pull required fields from deriv pack
                const auto Xv = xmmbar.X(ir, iz);
                const auto thXv = xmmbar.P_X(ir, iz);      // thorn X
                const auto edXv = xmmbar.EH_X(ir, iz);     // edthH X
                const auto thEdX = xmmbar.EH_P_X(ir, iz);   // edthH(thorn X) by commutation (this is what you used)

                // term1 = -1/2 * rhob * ( thorn(edthX) + (rho - rhob) edthX )
                const auto term1 = -half * rhob * (thEdX + (rho - rhob) * edXv);

                // term2 = -1/2 * rhob * (3 rhob + rho) * tau0 * (thorn X)
                const auto term2 = -half * rhob * ((three * rhob + rho) * tau0 * thXv);

                // term3 = +1/2 * 4 rho rhob^2 tau0 * X
                const auto term3 = half * (four * rho * rhob * rhob) * tau0 * Xv;

                const auto N_of_xmmbar = term1 + term2 + term3;

                Complex rhs = N_of_xmmbar.value();
                if (Tlm) rhs += (*Tlm)(ir, iz).value();

                NX_source.set_index(ir, iz, GHPScalar<Complex>(rhs, out_type.p, out_type.q));
            }
        }

        return NX_source;
    }

// -----------------------------------------------------------------------------
// x_nn source: T_ln + Re(U[x_mmbar]) + Re(V[x_nm])
// -----------------------------------------------------------------------------

    GHPSpectral source_xnn(
            const XDerivs &xmmbar,
            const XDerivs &xnm,
            const std::vector<Real> &r_grid,
            const std::vector<Real> &z_grid,
            const KerrMetricOutgoing &metric,
            const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords> &held,
            const ghp::GHPBackgroundFieldsVectorized &ghp,
            GHPType out_type,
            const GHPSpectral *Tln) {
        const size_t Nr = xmmbar.X.Nr();
        const size_t Nz = xmmbar.X.Nz();

        // Mode labels: assume match
        const auto modes = xmmbar.X.modes(); // adjust accessor if needed

        GHPSpectral source(
                Nr, Nz, modes,
                GHPScalar<Complex>(teuk::zeroC, out_type.p, out_type.q));

        const Real a = metric.a();

#pragma omp parallel for collapse(2) default(none) \
    shared(source, xmmbar, xnm, metric, held, ghp, Tln, a, Nr, Nz, teuk::zeroC, half, one, two, three, four, five, out_type)
        // Loop over grid points to assemble source using U and V formulas from GHZ reconstruction paper
        for (size_t iz = 0; iz < Nz; ++iz) {

            // Held-circle background on this z index (add what you need)
            const auto tau0 = held.tauH(iz);
            const auto taub0 = held.tauH_bar(iz);   // if you have it
            const auto psi0 = held.PsiH(iz);       // if you have it
            const auto rhoP0 = held.rhopH(iz);      // if you have it
            const auto Om0 = held.OmH(iz);        // if you have it

            for (size_t ir = 0; ir < Nr; ++ir) {

                // Background rho
                const auto rho_ghp = ghp.rho(ir, iz);
                const auto rhob_ghp = rho_ghp.conj();

                // optional external T_ln
                const Complex T = (Tln ? (*Tln)(ir, iz).value() : teuk::zeroC);

                // ---- Pull xmmbar derivs you might need for U ----
                const auto X = xmmbar.X(ir, iz);
                const auto PX = xmmbar.P_X(ir, iz);
                const auto EX = xmmbar.EH_X(ir, iz);
                const auto EbX = xmmbar.EbH_X(ir, iz);
                const auto EbEX = xmmbar.EbH_EH_X(ir, iz);
                const auto EPX = xmmbar.EH_P_X(ir, iz);
                const auto EbPX = xmmbar.EbH_P_X(ir, iz);
                const auto PpPX = xmmbar.PpHr_P_X(ir, iz);

                // ---- Pull xnm derivs you might need for V ----
                const auto Xn = xnm.X(ir, iz);
                const auto PXn = xnm.P_X(ir, iz);
                const auto EbXn = xnm.EbH_X(ir, iz);     // if you included it in xnm deriv pack
                const auto EbPXn = xnm.EbH_P_X(ir, iz);  // if you included it

                // -----------------------------------------------------------------
                // Implement U[xmmbar]
                // -----------------------------------------------------------------

                // --- Block A:  ρ \barρ [ - ð~'ð~ X  - \barτ°(ρ+ρ̄) ð~X + τ° ð~(P X) + τ°(ρ̄-ρ) ð~'X + τ° ð~'(P X) ]
                //
                // Here P ≡ thorn in outgoing/Kinnersley
                // So ð~(P X) is ð~(thorn X) which we already stored as thEdX.
                // And ð~'(P X) is ð~'(thorn X) which we stored as thEdbX.
                //
                const auto A =
                        (rho_ghp * rhob_ghp) * (
                                -1_r*(rho_ghp+rhob_ghp)*EX
                                - EbEX + tau0*EPX
                                + tau0*(rhob_ghp-rho_ghp)*EbX
                                + tau0*EbPX
                        );

                // --- Block B:  + \tilde P' P X
                const auto B = PpPX;
                // --- Block C:  + 1/2 [ ... ] P X  - 2 ρ \tilde P' X
                const auto bracket =
                        psi0*(rho_ghp*rho_ghp - rhob_ghp*rhob_ghp)
                        - two*rho_ghp*rhob_ghp*(rhob_ghp-rho_ghp)*tau0*taub0
                        - four*rho_ghp*rhob_ghp*rhoP0*(GHPScalar<Complex>(1.0, 0, 0) + rho_ghp*Om0);

                const auto C = half*bracket*PX - two*rho_ghp*PX; // here need to fix with  PH'
                // --- Block D: remaining purely multiplicative terms times X (no derivatives)
                const auto D =
                        (rho_ghp*rho_ghp) * psi0 * (two*rho_ghp-rhob_ghp)
                        + two * (rhob_ghp*rhob_ghp) * ( rhoP0-rho_ghp * (two*rho_ghp-rho_ghp) * tau0*taub0 )
                        - rho_ghp*rhob_ghp * (rho_ghp-three*rhob_ghp) * rhoP0*Om0;

                const auto UXmmb = (A+B+C+D);


                //
                // V[Xnm] = -rho * ( -3 rhob * EbX + (EbPX - tau0*(2 rhob + rho)*PX) + 2 rho^2 tau0 * X )
                //
                const auto VXnm =
                        -one*rho_ghp*(
                                -three*rhob_ghp*EbXn
                                + ( EbPXn-tau0 * (two*rhob_ghp+rho_ghp)*PXn )
                                + two * (rho_ghp*rho_ghp) * tau0*Xn
                        );

                const Real rhs_real = std::real(T) + std::real(UXmmb.value()) + std::real(VXnm.value());
                source.set_index(ir, iz, GHPScalar<Complex>(Complex(rhs_real, 0.0), out_type.p, out_type.q));
            }
        }

        return source;
    }
} //ghz::source