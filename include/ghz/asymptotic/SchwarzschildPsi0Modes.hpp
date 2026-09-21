#ifndef GHZ_ASYMPTOTIC_SCHWARZSCHILD_PSI0_MODES_HPP
#define GHZ_ASYMPTOTIC_SCHWARZSCHILD_PSI0_MODES_HPP

#include "ghz/core/GhzTypes.hpp"

#include <filesystem>
#include <map>
#include <utility>
#include <vector>

namespace ghz::asymptotic {

class SchwarzschildPsi0Modes {
public:
    static SchwarzschildPsi0Modes load_csv(const std::filesystem::path& path);

    [[nodiscard]] teuk::Complex coefficient(int ell, int m) const;
    [[nodiscard]] std::vector<teuk::Complex> summed_reduced_mode(
        int m, const std::vector<teuk::Real>& z) const;
    [[nodiscard]] std::vector<teuk::Complex> exact_reduced_seed(
        int m, teuk::Real omega, const std::vector<teuk::Real>& z) const;

    [[nodiscard]] int ell_min() const noexcept { return ell_min_; }
    [[nodiscard]] int ell_max() const noexcept { return ell_max_; }
    [[nodiscard]] teuk::Real mass() const noexcept { return mass_; }
    [[nodiscard]] teuk::Real orbital_radius() const noexcept { return orbital_radius_; }

    static teuk::Real pole_factor(int m, int spin_weight, teuk::Real z);
    static teuk::Real spin_weighted_spherical_harmonic(
        int spin_weight, int ell, int m, teuk::Real z);
    static teuk::Real reduced_spin_weighted_spherical_harmonic(
        int spin_weight, int ell, int m, teuk::Real z);
    static teuk::Real reduced_spin2_spherical_harmonic(int ell, int m, teuk::Real z);

private:
    int ell_min_ = 2;
    int ell_max_ = 0;
    teuk::Real mass_ = 1;
    teuk::Real orbital_radius_ = 10;
    std::map<std::pair<int, int>, teuk::Complex> coefficients_;
};

} // namespace ghz::asymptotic

#endif
