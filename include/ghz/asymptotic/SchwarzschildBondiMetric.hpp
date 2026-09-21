#ifndef GHZ_ASYMPTOTIC_SCHWARZSCHILD_BONDI_METRIC_HPP
#define GHZ_ASYMPTOTIC_SCHWARZSCHILD_BONDI_METRIC_HPP

#include "ghz/asymptotic/SchwarzschildPsi0Modes.hpp"

#include <array>

namespace ghz::asymptotic {

/** Coefficients of A_1 r + A_0 + A_{-1}/r + A_{-2}/r^2. */
struct RadialLaurent {
    std::array<teuk::Complex, 4> coefficient{};

    [[nodiscard]] teuk::Complex at(teuk::Real radius) const;
};

struct MetricComponentParts {
    RadialLaurent reconstructed;
    RadialLaurent lie_zeta;
    RadialLaurent adjusted;
};

struct SchwarzschildBondiMetricMode {
    MetricComponentParts nn;  // spin weight 0
    MetricComponentParts nm;  // spin weight +1
    MetricComponentParts mm;  // spin weight +2
};

/** Schwarzschild modal reference for 2 Re S-dagger(Phi^Adj) + L_zeta g.
 *
 * Follows the current BondiGauge_Adjusted.nb (mostly-minus signature),
 * including its plus sign and explicit factor of two. This is a closed-form
 * ell,m reference; the C++ m-mode Held collocation is checked separately.
 */
class SchwarzschildBondiMetric {
public:
    explicit SchwarzschildBondiMetric(const SchwarzschildPsi0Modes& psi0)
        : psi0_(psi0) {}

    [[nodiscard]] SchwarzschildBondiMetricMode mode(int ell, int m,
                                                    teuk::Real omega) const;
    [[nodiscard]] SchwarzschildBondiMetricMode summed_at_z(int m,
                                                           teuk::Real omega,
                                                           teuk::Real z) const;

private:
    const SchwarzschildPsi0Modes& psi0_;
};

} // namespace ghz::asymptotic

#endif
