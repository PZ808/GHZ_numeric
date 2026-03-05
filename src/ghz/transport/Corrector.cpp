// Created by Peter Zimmerman on 01.12.25.
//
#include "ghz/core/GhzTypes.hpp"
#include "ghz/transport/ODE.hpp"
#include "ghz/core/MathMacros.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/transport/Corrector.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"

#include <optional>
namespace ghz {

    using namespace teuk::literals;
    Complex I = teuk::I;

    /// Returned field lives on the same (r_grid, z_grid) and can be sampled in BuilderXnm.
    inline spectral::SpectralGHPVectorized build_source_for_xmmb(
            const spectral::SpectralGHPVectorized* Tll = nullptr
    ){
        // For the x_mmbar equation, the source is just T_{ll} (if you have it as a grid field)
        // If Tll is nullptr, this will just return a zero field on the same grid and with the same mode labels as Tll would have had.
        if (Tll) {
            return *Tll; // copy constructor
        } else {
            // Return a zero field with the same grid and mode labels as Tll would have had
            // You can customize the p,q weights as needed; here we just set them to zero for simplicity
            return spectral::SpectralGHPVectorized(Tll->r_grid(), Tll->z_grid(), Tll->mode_labels(),
                                                  GHPScalar<Complex>(teuk::zeroC, 0, 0));
        }

    }

/// Build RHS source for the x_nm equation from the already solved for x_mmbar data and external T_{lm}.
/// - sol_xmmbar[iz][ir] state layout: [ReX, ReX', ImX, ImX'] \n
/// - held provides tau^∘(z) (and anything else you may add later) \n
/// - if Tlm is provided, it is added pointwise to the returned field
///
/// Returned field lives on the same (r_grid, z_grid) and can be sampled in BuilderXnm.
    inline spectral::SpectralGHPVectorized build_source_for_xnm(
            const std::vector<std::vector<ode::StateVec>>& sol_xmmbar,
            const std::vector<Real>& r_grid,
            const std::vector<Real>& z_grid,
            const KerrMetricOutgoing& g_kerr,
            const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords>& held_coeffs,
            const ghp::GHPBackgroundFieldsVectorized& ghp_coeffs,
            const KinnersleyHeldOperators<OutgoingCoords>& kh_ops,
            int p_xmmbar, int q_xmmbar,          // weights for X_mmbar (metadata only)
            int p_out, int q_out,                // weights for output source (metadata only)
            int m, int kr, int kz,               // mode labels stored in GHPSpectral
            const spectral::SpectralGHPVectorized* Tlm = nullptr
    ) {
        const size_t Nr = r_grid.size();
        const size_t Nz = z_grid.size();

        // Allocate containers for X and thornX
        // to be filled with the data from the ODE state results
        // (note that thornX = dX/dr in outgoing chart and Kinnerley tetrad)
        GHPSpectral X      (Nr, Nz, {m, kr, kz},
                            GHPScalar<Complex>(teuk::zeroC, p_xmmbar,   q_xmmbar));
        GHPSpectral thornX (Nr, Nz, {m, kr, kz},
                            GHPScalar<Complex>(teuk::zeroC, p_xmmbar+1, q_xmmbar+1));

        // parallelize the loop to populate X and thornX from the ODE solution for x_mmbar
#pragma omp parallel for collapse(2) default(none) shared(Nr, Nz, X, thornX, sol_xmmbar, r_grid, z_grid, p_xmmbar, q_xmmbar)
        for (size_t iz = 0; iz < Nz; ++iz) {
            for (size_t ir = 0; ir < Nr; ++ir) {
                // populate X and thornX from the ODE solution for x_mmbar
                const auto& y = sol_xmmbar[iz][ir]; // y = {ReX, ReX', ImX, ImX'} at (r,z)
                const Complex Xval (y[0], y[2]);
                const Complex DrX  (y[1], y[3]);     // <-- stored radial derivative from mmbar ODE

                X.set_index(ir, iz, GHPScalar<Complex>(Xval, p_xmmbar, q_xmmbar));
                thornX.set_index(ir, iz, GHPScalar<Complex>(DrX,  p_xmmbar+1, q_xmmbar+1));
            }
        }

        // Compute edthX = \tilde{edth} X and thornEdthX = \thorn(\tilde{edth}X)
        //         using commutation: thornEdthX = \tilde{edth}(thornX) (they commute for type-D backgrounds, see Held 1974)
        //
        GHPSpectral edthX       (Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p_xmmbar,   q_xmmbar-2));
        GHPSpectral thornEdthX  (Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p_xmmbar+1, q_xmmbar-1)); // (p+1,q+1)->(p+1,q-1)

#pragma omp parallel default(none) shared(g_kerr, X, kh_ops, thornX, edthX, thornEdthX, r_grid, z_grid, Nr, Nz)
        {
            spectral::SpectralDiffer differ(static_cast<int>(Nz), static_cast<int>(Nr));

#pragma omp for
            for (size_t ir = 0; ir < Nr; ++ir) {
                auto Xslice         = X.slice_R(ir);
                auto thornXslice    = thornX.slice_R(ir);

                auto edthXslice     = edthX.slice_R(ir);
                auto edthThXslice   = thornEdthX.slice_R(ir);

                kh_ops.edthH_inplace_RSliceV(Xslice,      edthXslice,   g_kerr.a());
                kh_ops.edthH_inplace_RSliceV(thornXslice, edthThXslice, g_kerr.a());
            }
        }

        GHPSpectral NX_source(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(
                teuk::zeroC, p_out, q_out));

#pragma omp parallel for collapse(2) default(none) \
        shared(Nr, Nz, NX_source, X, thornX, edthX, thornEdthX, r_grid, z_grid,I,\
        g_kerr, held_coeffs, ghp_coeffs, Tlm)
        for (int iz = 0; iz < Nz; ++iz) {
            // tau^∘ depends only on z (Held)
            auto tau0 = held_coeffs.tauH(iz);

            for (int ir = 0; ir < Nr; ++ir) {
                const Real r = r_grid[ir];
                const Real z = z_grid[iz];

                auto rho_ghp = ghp_coeffs.rho(ir, iz);
                auto rhobar_ghp = rho_ghp.conj();

                auto rho  = -Real(1)/(r-I*(0.0, g_kerr.a() * z));
                auto rhob = std::conj(rho);

                auto Xv    = X(ir, iz);
                auto thXv  = thornX(ir, iz);
                auto edXv  = edthX(ir, iz);
                auto thEdX = thornEdthX(ir, iz);

                // term1 = rhob * (thorn + rho - rhob) (edthX)
                auto term1 = -Real(0.5) * rhob * ( thEdX + (rho - rhob) * edXv );

                // term2 = - rhob * (3 rhob + rho) * tau0 * (thorn X)
                auto term2 = - Real(0.5)*rhob * ( (Real(3) * rhob + rho) * tau0 * thXv );

                // term3 = + 4 rho rhob^2 tau0 * X
                auto term3 = Real(0.5)*(Real(4) * rho * rhob * rhob) * tau0 * Xv;

                auto N_of_xmmbar = term1+term2+term3; //  (term1 + term2 + term3);

                Complex rhs = N_of_xmmbar.value();
                if (Tlm) rhs += (*Tlm)(ir, iz).value();  // add T_lm if you have it as a grid field

                NX_source.set_index(ir, iz, GHPScalar<Complex>(rhs, NX_source.p(), NX_source.q()));
            }
        }

        return NX_source;
    }


    inline spectral::SpectralGHPVectorized build_source_for_xnn(
            const std::vector<std::vector<ode::StateVec>>& sol_xmmbar,
            const std::vector<std::vector<ode::StateVec>>& sol_xnm,
            const std::vector<Real>& r_grid,
            const std::vector<Real>& z_grid,
            const KerrMetricOutgoing& g_kerr,
            const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords>& held_coeffs,
            const ghp::GHPBackgroundFieldsVectorized& ghp_coeffs,
            const KinnersleyHeldOperators<OutgoingCoords>& kh_ops,
            int p_xmmbar, int q_xmmbar,
            int p_xnm,    int q_xnm,
            int p_out,    int q_out,
            int m, int kr, int kz,
            const spectral::SpectralGHPVectorized* Tln /*= nullptr*/)
    {

        const size_t Nr = r_grid.size();
        const size_t Nz = z_grid.size();
        const Real a = g_kerr.a();

        // -----------------------------
        // 1) Reconstruct lower-level fields + stored thorn from previous ODE states
        // -----------------------------
        GHPSpectral Xmmb(Nr, Nz, {m, kr, kz},
                         GHPScalar<Complex>(teuk::zeroC, p_xmmbar,     q_xmmbar));
        GHPSpectral P_Xmmb(Nr, Nz, {m, kr, kz},
                           GHPScalar<Complex>(teuk::zeroC, p_xmmbar + 1, q_xmmbar + 1));

        GHPSpectral Xnm(Nr, Nz, {m, kr, kz},
                        GHPScalar<Complex>(teuk::zeroC, p_xnm,        q_xnm));
        GHPSpectral P_Xnm(Nr, Nz, {m, kr, kz},
                          GHPScalar<Complex>(teuk::zeroC, p_xnm + 1,    q_xnm + 1));

#pragma omp parallel for collapse(2) default(none) \
        shared(Nr, Nz, sol_xmmbar, sol_xnm, Xmmb, P_Xmmb, Xnm, P_Xnm, p_xmmbar, q_xmmbar, p_xnm, q_xnm)
        for (size_t iz = 0; iz < Nz; ++iz) {
            for (size_t ir = 0; ir < Nr; ++ir) {
                // xmmbar: y = {ReX, ReX', ImX, ImX'}
                {
                    const auto& y = sol_xmmbar[iz][ir];
                    const Complex Xval (y[0], y[2]);
                    const Complex DrX  (y[1], y[3]);

                    Xmmb.set_index(ir, iz, GHPScalar<Complex>(Xval, p_xmmbar, q_xmmbar));
                    P_Xmmb.set_index(ir, iz, GHPScalar<Complex>(DrX, p_xmmbar + 1, q_xmmbar + 1));
                }
                // xnm: same packing assumption
                {
                    const auto& y = sol_xnm[iz][ir];
                    const Complex Xval (y[0], y[2]);
                    const Complex DrX  (y[1], y[3]);

                    Xnm.set_index(ir, iz, GHPScalar<Complex>(Xval, p_xnm, q_xnm));
                    P_Xnm.set_index(ir, iz, GHPScalar<Complex>(DrX, p_xnm + 1, q_xnm + 1));
                }
            }
        } // populate Xmmbar, thornXmmbar, Xnm, thornXnm from ODE states for x_mmbar and x_nm

        // -----------------------------
        // 2) Held angular derivatives on each r-slice
        //    You need:  ð~,  ð~'  and also ð~'ð~ acting on xmmbar,
        //    and for the V-operator you need ð~' acting on xnm plus P acting on xnm.
        // -----------------------------

        // xmmbar derivatives
        GHPSpectral EH_Xmmb     (Nr, Nz, {m, kr, kz},
                                 GHPScalar<Complex>(teuk::zeroC, p_xmmbar,     q_xmmbar - 2));
        GHPSpectral EbH_Xmmb   (Nr, Nz, {m, kr, kz},
                                GHPScalar<Complex>(teuk::zeroC, p_xmmbar - 2,   q_xmmbar));
        GHPSpectral EH_EbH_Xmmb (Nr, Nz, {m, kr, kz},
                                 GHPScalar<Complex>(teuk::zeroC, p_xmmbar - 2, q_xmmbar - 2));
        GHPSpectral EbH_EH_Xmmb (Nr, Nz, {m, kr, kz},
                                 GHPScalar<Complex>(teuk::zeroC, p_xmmbar - 2, q_xmmbar - 2));

        GHPSpectral PpH_P_Xmmb   (Nr, Nz, {m, kr, kz},
                                  GHPScalar<Complex>(teuk::zeroC, p_xmmbar, q_xmmbar));
        GHPSpectral EH_P_Xmmbar  (Nr, Nz, {m, kr, kz},
                                  GHPScalar<Complex>(teuk::zeroC, p_xmmbar + 1, q_xmmbar - 1));
        GHPSpectral EbH_P_Xmmbar (Nr, Nz, {m, kr, kz},
                                  GHPScalar<Complex>(teuk::zeroC, p_xmmbar - 1, q_xmmbar + 1));

        // xnm derivatives needed by V
        GHPSpectral EbH_Xnm   (Nr, Nz, {m, kr, kz},
                               GHPScalar<Complex>(teuk::zeroC, p_xnm - 2,   q_xnm));
        GHPSpectral EbH_P_Xnm (Nr, Nz, {m, kr, kz},
                               GHPScalar<Complex>(teuk::zeroC, p_xnm - 1,   q_xnm + 1));

#pragma omp parallel default(none) \
        shared(Nr, Nz, a, kh_ops, Xmmb,P_Xmmb, EbH_Xmmb, EH_Xmmb, EH_EbH_Xmmb, \
        EbH_EH_Xmmb, PpH_P_Xmmb, EH_P_Xmmbar, \
               EbH_P_Xmmbar, Xnm, P_Xnm, EbH_Xnm, EbH_P_Xnm)
        {

#pragma omp for
            for (size_t ir=0; ir<Nr; ++ir) {
                //
                // --- Xmmbar ---
                //
                auto Xmmb_slice       = Xmmb.slice_R(ir);
                auto P_Xmmb_slice     = P_Xmmb.slice_R(ir);
                auto PpH_P_Xmmb_slice = PpH_P_Xmmb.slice_R(ir);
                auto EH_Xmmb_slice    = EH_Xmmb.slice_R(ir);
                auto EbH_Xmmbar_slice = EbH_Xmmb.slice_R(ir);
                auto EbH_EH_Xmmb_slice = EbH_EH_Xmmb.slice_R(ir);
                auto EH_EbH_Xmmb_slice = EH_EbH_Xmmb.slice_R(ir);
                auto EH_P_Xmmb_slice   = EH_P_Xmmbar.slice_R(ir);
                auto EbH_P_Xmmb_slice  = EbH_P_Xmmbar.slice_R(ir);

                kh_ops.thornPH_inplace_RSliceV(P_Xmmb_slice, PpH_P_Xmmb_slice);
                kh_ops.edthH_inplace_RSliceV(Xmmb_slice,   EH_Xmmb_slice);
                kh_ops.edthBarH_inplace_RSliceV(Xmmb_slice, EbH_Xmmbar_slice);
                kh_ops.edthBarH_inplace_RSliceV(EH_Xmmb_slice, EbH_EH_Xmmb_slice);
                kh_ops.edthH_inplace_RSliceV(P_Xmmb_slice,   EH_P_Xmmb_slice);
                kh_ops.edthBarH_inplace_RSliceV(P_Xmmb_slice, EbH_P_Xmmb_slice);
                //
                // --- Xnm ---
                //
                auto Xnm_slice   = Xnm.slice_R(ir);
                auto P_Xnm_slice = P_Xnm.slice_R(ir);

                auto EbH_Xnm_slice   = EbH_Xnm.slice_R(ir);
                auto EbH_P_Xnm_slice = EbH_P_Xnm.slice_R(ir);

                // Need ð~' X_nm and ð~'(thorn X_nm)
                kh_ops.edthBarH_inplace_RSliceV(Xnm_slice,   EbH_Xnm_slice);
                kh_ops.edthBarH_inplace_RSliceV(P_Xnm_slice, EbH_P_Xnm_slice);
            }
        }  // loop over r-slices to compute Held derivatives needed for source
        // -----------------------------
        // 3) Assemble RHS pointwise on the Nr x Nz grid using the Held formulas for U and V from the paper,
        // plus T_{ln} if you have it as a grid function.
        //     RHS = T_ln + Re[ U * x_mmbar ] + Re[ V * x_nm ].
        // -----------------------------
        GHPSpectral source(Nr, Nz, {m, kr, kz},
                           GHPScalar<Complex>(teuk::zeroC, p_out, q_out)); // container for the final source for x_nn

        auto get_Tln = [&](size_t ir, size_t iz) -> Complex {
            if (!Tln) return teuk::zeroC;
            return (*Tln)(ir, iz).value();
        };


        // Apply U to xmmbar at a point (ir,iz)
        auto apply_U_to_xmmbar = [&](size_t ir, size_t iz) -> Complex {
            // field values
            auto X    = Xmmb(ir, iz);
            auto PX   = P_Xmmb(ir, iz);
            auto EX   = EH_Xmmb(ir, iz);           // ð~ X
            auto EbX  = EbH_Xmmb(ir, iz);           // ð~ X
            auto EbEX = EbH_EH_Xmmb(ir, iz);    // ð~'ð~ X
            auto PpPX = PpH_P_Xmmb(ir,iz); // thorn~'(thorn X)
            auto EPX  = EH_P_Xmmbar(ir, iz);      // thorn(ð~X)=ð~(thorn X)
            auto EbPX = EbH_P_Xmmbar(ir, iz);   // thorn(ð~'X)=ð~'(thorn X)
            // background at this z-index (vectorized in z)
            // NOTE: replace these with *your actual accessors*.

            auto tau0    = held_coeffs.tauH(iz);        // τ^°
            auto taub0   = held_coeffs.tauH_bar(iz);     // \bar τ^°
            auto psi0    = held_coeffs.PsiH(iz);        // Ψ^°
            auto rhoP0   = held_coeffs.rhopH(iz);   // ρ'^°
            auto Omega0  = held_coeffs.OmH(iz);      // Ω^°

            // geometric rho(r,z)
            // NOTE: use your existing rho(r,z) helper if you have it; I’m writing it explicitly.
            // rho = -1/(r - i a z)
            // rhob = -1/(r + i a z)
            // (check sign conventions vs your codebase)
            const Real r = r_grid[ir];
            const Real z = z_grid[iz];
            const Complex denom (r, -a*z);
            //const GHPScalar rho   = GHPScalar(-Complex(1.0,0.0)/denom,1,1);
            //const GHPScalar rhob  = GHPScalar(-Complex(1.0,0.0)/conj(denom),1,1);

            auto rho_ghp = ghp_coeffs.rho(ir, iz);
            auto rhob_ghp= rho_ghp.conj();

            // Now build U[X] from the paper formula.
            // I’m writing it in “term blocks” so it’s easy to line up with the PDF.

            // --- Block A:  ρ \barρ [ - ð~'ð~ X  - \barτ°(ρ+ρ̄) ð~X + τ° ð~(P X) + τ°(ρ̄-ρ) ð~'X + τ° ð~'(P X) ]
            //
            // Here P ≡ thorn in your outgoing/Kinnersley setup.
            // So ð~(P X) is ð~(thorn X) which we already stored as thEdX.
            // And ð~'(P X) is ð~'(thorn X) which we stored as thEdbX.
            //
            auto A =
                    (rho_ghp*rhob_ghp) * (
                            -1_r*(rho_ghp+rhob_ghp)*EX-EbEX + tau0*EPX
                            + tau0*(rhob_ghp-rho_ghp)*EbX
                            + tau0*EbPX
                    );

            // --- Block B:  + \tilde P' P X
            auto B = PpH_P_Xmmb(ir, iz);

            // --- Block C:  + 1/2 [ ... ] P X  - 2 ρ \tilde P' X

            auto bracket =
                    psi0*(rho_ghp*rho_ghp - rhob_ghp*rhob_ghp)
                    - two*rho_ghp*rhob_ghp*(rhob_ghp-rho_ghp)*tau0*taub0
                    - four*rho_ghp*rhob_ghp*rhoP0*(GHPScalar<Complex>(1.0,0,0) + rho_ghp*Omega0);

            auto C = half*bracket*PX -two*rho_ghp*PX; // here need to fix with  PH'

            // --- Block D: remaining purely multiplicative terms times X (no derivatives)

            auto D =
                    (rho_ghp*rho_ghp)*psi0*(two*rho_ghp-rhob_ghp)
                    + two*(rhob_ghp*rhob_ghp)*( rhoP0 - rho_ghp*(two*rho_ghp-rho_ghp)*tau0*taub0 )
                    - rho_ghp*rhob_ghp*(rho_ghp-three*rhob_ghp)*rhoP0*Omega0;

            return (A + B + C + D).value();
        }; // apply U to x_mmbar at a point (ir,iz)

        // Apply V to xnm at a point (ir,iz)
        auto apply_V_to_xnm = [&](size_t ir, size_t iz) -> Complex {
            auto X   = Xnm(ir, iz);
            auto PX = P_Xnm(ir, iz);

            auto EbX   = EbH_Xnm(ir, iz);      // ð~' X
            auto EbPX = EbH_P_Xnm(ir, iz); // ð~'(thorn X)

            auto tau0 = held_coeffs.tauH(iz);

            const Real r = r_grid[ir];
            const Real z = z_grid[iz];
            auto rho_ghp = ghp_coeffs.rho(ir, iz);
            auto rhob_ghp= rho_ghp.conj();

            //
            // operator V
            //   V = -ρ{ -3 ρ̄ ð~' + [ ð~' - τ°(2 ρ̄ + ρ)] P + 2 ρ^2 τ° }
            // Acting on X:
            //   V[X] = -ρ( -3 ρ̄ ð~'X + (ð~'(P X) - τ°(2ρ̄+ρ) P X) + 2 ρ^2 τ° X )
            //
            auto VX =
                    -one*rho_ghp * (
                            -three*rhob_ghp*EbX
                            + (EbPX - tau0*(two*rhob_ghp + rho_ghp)*PX)
                            + two*(rho_ghp*rho_ghp)*tau0*X
                    );

            return VX.value();
        };

#pragma omp parallel for collapse(2) default(none) \
        shared(Nr, Nz, source, Tln, get_Tln, apply_U_to_xmmbar, apply_V_to_xnm)
        for (size_t iz = 0; iz < Nz; ++iz) {
            for (size_t ir = 0; ir < Nr; ++ir) {
                const Complex T = get_Tln(ir, iz);

                const Complex Ux = apply_U_to_xmmbar(ir, iz);
                const Complex Vx = apply_V_to_xnm(ir, iz);

                // RHS uses Re[...] pieces:
                const Real rhs_real = std::real(T) + std::real(Ux) + std::real(Vx);

                source.set_index(ir, iz, GHPScalar<Complex>(Complex(rhs_real, 0.0), p_out, q_out));
            }
        }

        return source;
    }

} // namespace ghz