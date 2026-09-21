#ifndef GHZ_ASYMPTOTIC_BONDI_HELD_SEED_SOLVE_HPP
#define GHZ_ASYMPTOTIC_BONDI_HELD_SEED_SOLVE_HPP

#include "ghz/core/GhzTypes.hpp"
#include "ghz/geom/Coords.hpp"
#include "ghz/spectral/KinnersleySpectralHeldOperators.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"

#include <vector>

namespace ghz::asymptotic {

struct BondiHeldSeedSolution {
    std::vector<teuk::Complex> reduced_seed;
    teuk::Real relative_residual = 0;
    teuk::Real equilibrated_condition_number = 0;
};

class BondiHeldSeedSolver {
public:
    BondiHeldSeedSolver(
        const spectral::SpectralDiffer& differ,
        const KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>& held_operators,
        teuk::Real mass);

    [[nodiscard]] std::vector<teuk::Complex> apply(
        const std::vector<teuk::Complex>& reduced_seed, int m, teuk::Real omega) const;

    [[nodiscard]] BondiHeldSeedSolution solve(
        const std::vector<teuk::Complex>& reduced_psi0, int m, teuk::Real omega) const;

private:
    const spectral::SpectralDiffer& differ_;
    const KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>& held_operators_;
    teuk::Real mass_;
};

} // namespace ghz::asymptotic

#endif
