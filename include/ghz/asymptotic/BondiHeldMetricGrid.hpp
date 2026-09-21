#ifndef GHZ_ASYMPTOTIC_BONDI_HELD_METRIC_GRID_HPP
#define GHZ_ASYMPTOTIC_BONDI_HELD_METRIC_GRID_HPP

#include "ghz/core/GhzTypes.hpp"
#include "ghz/geom/Coords.hpp"
#include "ghz/spectral/KinnersleySpectralHeldOperators.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"

#include <array>
#include <vector>

namespace ghz::asymptotic {

struct MetricGridLaurent {
    // Coefficients of r, 1, 1/r, and 1/r^2 at each z node.
    std::array<std::vector<teuk::Complex>, 4> coefficient;
};

struct MetricGridParts {
    MetricGridLaurent reconstructed;
    MetricGridLaurent lie_zeta;
    MetricGridLaurent adjusted;
};

struct BondiHeldMetricGrid {
    MetricGridParts nn; // spin 0
    MetricGridParts nm; // spin +1
    MetricGridParts mm; // spin +2
};

/** Schwarzschild m-mode Held-grid reconstruction of the adjustment metric.
 *
 * Inputs are reduced spin +2 f_m and reduced spin -2 conjugate fbar_m on the
 * same LGL z grid. All angular operations use the pole-factorized spectral
 * Held operators, without an ell projection.
 */
class BondiHeldMetricReconstruction {
public:
    BondiHeldMetricReconstruction(
        const spectral::SpectralDiffer& differ,
        const KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>& operators,
        teuk::Real mass);

    [[nodiscard]] BondiHeldMetricGrid reconstruct(
        const std::vector<teuk::Complex>& reduced_f,
        const std::vector<teuk::Complex>& reduced_fbar,
        int m, teuk::Real omega) const;

private:
    [[nodiscard]] std::vector<teuk::Complex> up(
        const std::vector<teuk::Complex>& value, int spin,
        int m, teuk::Real omega) const;
    [[nodiscard]] std::vector<teuk::Complex> down(
        const std::vector<teuk::Complex>& value, int spin,
        int m, teuk::Real omega) const;

    const spectral::SpectralDiffer& differ_;
    const KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>& operators_;
    teuk::Real mass_;
};

} // namespace ghz::asymptotic

#endif
