//
// TeukolskySource.hpp
//
// Assembly of the spin +2 Teukolsky effective source from tetrad components.
//

#ifndef GHZ_SOURCE_TEUKOLSKYSOURCE_HPP
#define GHZ_SOURCE_TEUKOLSKYSOURCE_HPP
#pragma once

#include "ghz/geom/KinnersleyTetrad.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"

#include <functional>
#include <stdexcept>
#include <vector>

namespace ghz::source {

    using teuk::Real;
    using teuk::Complex;

    using TeukolskyField = spectral::SpectralGHPVectorized;

    struct TeukolskySourceOperators {
        using UnaryOperator =
                std::function<void(const TeukolskyField&, TeukolskyField&)>;

        UnaryOperator thorn;
        UnaryOperator eth;
    };

    struct TeukolskySourceCoefficients {
        TeukolskyField rho;
        TeukolskyField rho_bar;
        TeukolskyField tau;
        TeukolskyField tau_prime_bar;
    };

    struct TeukolskySourceComponents {
        TeukolskyField Tll;
        TeukolskyField Tlm_sym;
        TeukolskyField Tmm;
    };

    // Implements
    //
    // S_0 T = 1/2 (eth - taupbar - 4 tau)
    //              [ (thorn - 2 rhobar) T_(lm)
    //                - (eth - taupbar) T_ll ]
    //       + 1/2 (thorn - 4 rho - rhobar)
    //              [ (eth - 2 taupbar) T_(lm)
    //                - (thorn - rhobar) T_mm ].
    //
    // The coefficient fields are pointwise background fields on the same
    // collocation grid. Their mode metadata is not used; the returned field
    // keeps the mode metadata and omega of the source components.
    [[nodiscard]] TeukolskyField apply_s0_teukolsky_source(
            const TeukolskySourceComponents& components,
            const TeukolskySourceCoefficients& coeffs,
            const TeukolskySourceOperators& ops);

    // Samples the Kinnersley tetrad's GHP spin coefficients on a tensor-product
    // coordinate grid and extracts the coefficients required by S_0 T.
    //
    // For the current BLCoords and OutgoingCoords implementations, x1_nodes are
    // r values and z_nodes are cos(theta). For OutgoingCoordsCompact, x1_nodes
    // are compact sigma values.
    template <typename CoordT>
    [[nodiscard]] TeukolskySourceCoefficients
    make_teukolsky_source_coefficients_from_tetrad(
            KinnersleyTetrad<CoordT>& tetrad,
            const std::vector<Real>& x1_nodes,
            const std::vector<Real>& z_nodes,
            spectral::SpectralFieldVectorized<ghp::GHPScalar<Complex>>::Modes modes =
                    {0, 0, 0}) {
        if (x1_nodes.empty() || z_nodes.empty()) {
            throw std::runtime_error(
                    "make_teukolsky_source_coefficients_from_tetrad: empty grid");
        }

        auto make_field = [&](int p, int q) {
            return TeukolskyField(x1_nodes.size(),
                                  z_nodes.size(),
                                  modes,
                                  ghp::GHPScalar<Complex>(teuk::zeroC, p, q),
                                  p,
                                  q);
        };

        TeukolskySourceCoefficients coeffs{
                make_field(1, 1),
                make_field(1, 1),
                make_field(1, -1),
                make_field(1, -1)};

        for (size_t ir = 0; ir < x1_nodes.size(); ++ir) {
            for (size_t iz = 0; iz < z_nodes.size(); ++iz) {
                CoordT X(Real(0), x1_nodes[ir], z_nodes[iz], Real(0));
                tetrad.build_tetrad_at(X);
                const auto scalars = tetrad.get_scalars_at(X).ghp_scalars;

                coeffs.rho(ir, iz) = scalars.rho;
                coeffs.rho_bar(ir, iz) = scalars.rho_bar;
                coeffs.tau(ir, iz) = scalars.tau;
                coeffs.tau_prime_bar(ir, iz) = scalars.taup_bar;
            }
        }

        return coeffs;
    }

    [[nodiscard]] inline TeukolskySourceCoefficients
    make_bl_teukolsky_source_coefficients_from_theta(
            KinnersleyTetrad<BLCoords>& tetrad,
            const std::vector<Real>& r_nodes,
            const std::vector<Real>& theta_nodes,
            spectral::SpectralFieldVectorized<ghp::GHPScalar<Complex>>::Modes modes =
                    {0, 0, 0}) {
        std::vector<Real> z_nodes(theta_nodes.size());
        for (size_t i = 0; i < theta_nodes.size(); ++i) {
            z_nodes[i] = teuk::Cos(theta_nodes[i]);
        }

        return make_teukolsky_source_coefficients_from_tetrad(
                tetrad,
                r_nodes,
                z_nodes,
                modes);
    }

} // namespace ghz::source

#endif // GHZ_SOURCE_TEUKOLSKYSOURCE_HPP
