#include "ghz/asymptotic/SchwarzschildPsi0Modes.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ghz::asymptotic {
namespace {

teuk::Real jacobi_polynomial(int n, int alpha, int beta, teuk::Real z) {
    if (n == 0) return teuk::Real(1);
    teuk::Real pm2 = 1;
    teuk::Real pm1 = (teuk::Real(alpha - beta) + teuk::Real(alpha + beta + 2) * z) / 2;
    if (n == 1) return pm1;
    for (int k = 2; k <= n; ++k) {
        const teuk::Real kk = k;
        const teuk::Real ab = alpha + beta;
        const teuk::Real lhs = 2 * kk * (kk + ab) * (2 * kk + ab - 2);
        const teuk::Real a = (2 * kk + ab - 1) *
            ((2 * kk + ab) * (2 * kk + ab - 2) * z + alpha * alpha - beta * beta);
        const teuk::Real b = 2 * (kk + alpha - 1) * (kk + beta - 1) * (2 * kk + ab);
        const teuk::Real p = (a * pm1 - b * pm2) / lhs;
        pm2 = pm1;
        pm1 = p;
    }
    return pm1;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    return fields;
}

teuk::Real parse_real(const std::string& text) {
    return static_cast<teuk::Real>(std::stold(text));
}

int parity_sign(int value) { return (std::abs(value) % 2 == 0) ? 1 : -1; }

} // namespace

SchwarzschildPsi0Modes SchwarzschildPsi0Modes::load_csv(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("cannot open psi0 coefficient table: " + path.string());

    SchwarzschildPsi0Modes result;
    bool mostly_minus = false;
    std::string line;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') {
            const auto equals = line.find('=');
            if (equals == std::string::npos) continue;
            const std::string key = line.substr(2, equals - 2);
            const std::string value = line.substr(equals + 1);
            if (key == "signature") mostly_minus = (value == "mostly-minus");
            else if (key == "mass") result.mass_ = parse_real(value);
            else if (key == "orbital_radius") result.orbital_radius_ = parse_real(value);
            else if (key == "ell_min") result.ell_min_ = std::stoi(value);
            else if (key == "ell_max") result.ell_max_ = std::stoi(value);
            continue;
        }
        if (line.rfind("ell,", 0) == 0) continue;
        const auto fields = split_csv(line);
        if (fields.size() != 4) throw std::runtime_error("malformed psi0 coefficient row");
        const int ell = std::stoi(fields[0]);
        const int m = std::stoi(fields[1]);
        result.coefficients_[{ell, m}] = {parse_real(fields[2]), parse_real(fields[3])};
    }
    if (!mostly_minus)
        throw std::runtime_error("psi0 coefficient table is not marked mostly-minus");
    const std::size_t expected = static_cast<std::size_t>(
        (result.ell_max_ + 1) * (result.ell_max_ + 1) - result.ell_min_ * result.ell_min_);
    if (result.coefficients_.size() != expected)
        throw std::runtime_error("psi0 coefficient table has incomplete (ell,m) coverage");
    return result;
}

teuk::Complex SchwarzschildPsi0Modes::coefficient(int ell, int m) const {
    const auto it = coefficients_.find({ell, m});
    if (it == coefficients_.end()) throw std::out_of_range("psi0 (ell,m) coefficient is absent");
    return it->second;
}

teuk::Real SchwarzschildPsi0Modes::pole_factor(int m, int spin_weight, teuk::Real z) {
    return std::pow(teuk::Real(1) - z, teuk::Real(std::abs(m + spin_weight)) / 2) *
           std::pow(teuk::Real(1) + z, teuk::Real(std::abs(m - spin_weight)) / 2);
}

teuk::Real SchwarzschildPsi0Modes::reduced_spin_weighted_spherical_harmonic(
    int spin_weight, int ell, int m, teuk::Real z) {
    if (std::abs(spin_weight) > ell || std::abs(m) > ell)
        throw std::invalid_argument("spin-weighted harmonic requires |s|, |m| <= ell");
    if (z < teuk::Real(-1) || z > teuk::Real(1))
        throw std::invalid_argument("spin-weighted harmonic requires -1 <= z <= 1");

    // _sY_lm(theta,0) = (-1)^s sqrt((2l+1)/(4pi)) d^l_{m,-s}(theta).
    // Remove the known pole factor from the Wigner-d representation.  The
    // remaining Jacobi polynomial is regular at both LGL endpoints.
    const int m_prime = m;
    const int m_second = -spin_weight;
    const int k_min = std::max(0, m_second - m_prime);
    const int k_max = std::min(ell + m_second, ell - m_prime);
    const int north_power = std::abs(m - spin_weight);
    const int south_power = std::abs(m + spin_weight);

    const teuk::Real log_prefactor = teuk::Real(0.5) * (
        std::lgamma(teuk::Real(ell + m_second + 1))
        + std::lgamma(teuk::Real(ell - m_second + 1))
        + std::lgamma(teuk::Real(ell + m_prime + 1))
        + std::lgamma(teuk::Real(ell - m_prime + 1)));
    const teuk::Real normalization = parity_sign(spin_weight)
        * std::sqrt(teuk::Real(2 * ell + 1) /
                    (teuk::Real(4) * std::numbers::pi_v<teuk::Real>));
    const teuk::Real two_to_minus_ell =
        std::pow(teuk::Real(2), teuk::Real(-ell));

    // Fix the normalization at the north pole, where the pole-factorized
    // Wigner sum contains a single nonzero term.  Evaluate the angular
    // dependence with the stable Jacobi recurrence instead of the strongly
    // cancelling Wigner sum.
    teuk::Real north_pole_value = 0;
    for (int k = k_min; k <= k_max; ++k) {
        const int cos_power = 2 * ell + m_second - m_prime - 2 * k;
        const int sin_power = m_prime - m_second + 2 * k;
        const int reduced_north_power = (cos_power - north_power) / 2;
        const int reduced_south_power = (sin_power - south_power) / 2;
        const teuk::Real log_denominator =
            std::lgamma(teuk::Real(ell + m_second - k + 1))
            + std::lgamma(teuk::Real(k + 1))
            + std::lgamma(teuk::Real(m_prime - m_second + k + 1))
            + std::lgamma(teuk::Real(ell - m_prime - k + 1));
        const teuk::Real coefficient = parity_sign(m_prime - m_second + k)
            * std::exp(log_prefactor - log_denominator) * two_to_minus_ell;
        if (reduced_south_power == 0) {
            north_pole_value += coefficient
                * std::pow(teuk::Real(2), reduced_north_power);
        }
    }
    const int degree = ell - std::max(std::abs(m), std::abs(spin_weight));
    const teuk::Real jacobi_at_north =
        std::exp(std::lgamma(teuk::Real(degree + south_power + 1))
                 - std::lgamma(teuk::Real(degree + 1))
                 - std::lgamma(teuk::Real(south_power + 1)));
    return normalization * north_pole_value
        * jacobi_polynomial(degree, south_power, north_power, z)
        / jacobi_at_north;
}

teuk::Real SchwarzschildPsi0Modes::spin_weighted_spherical_harmonic(
    int spin_weight, int ell, int m, teuk::Real z) {
    return pole_factor(m, spin_weight, z)
        * reduced_spin_weighted_spherical_harmonic(spin_weight, ell, m, z);
}

teuk::Real SchwarzschildPsi0Modes::reduced_spin2_spherical_harmonic(
    int ell, int m, teuk::Real z) {
    if (m < 2 || ell < m)
        throw std::invalid_argument("reduced spin-2 harmonic currently requires ell >= m >= 2");
    const teuk::Real log_ratio = std::lgamma(teuk::Real(ell + m + 1)) +
        std::lgamma(teuk::Real(ell - m + 1)) - std::lgamma(teuk::Real(ell - 1)) -
        std::lgamma(teuk::Real(ell + 3));
    const teuk::Real normalization = parity_sign(m - 2) *
        std::sqrt(teuk::Real(2 * ell + 1) / (4 * std::numbers::pi_v<teuk::Real>)) *
        std::exp(log_ratio / 2) / std::pow(teuk::Real(2), m);
    const teuk::Real result =
        normalization * jacobi_polynomial(ell - m, m + 2, m - 2, z);
    return result;
}

std::vector<teuk::Complex> SchwarzschildPsi0Modes::summed_reduced_mode(
    int m, const std::vector<teuk::Real>& z) const {
    std::vector<teuk::Complex> values(z.size(), teuk::zeroC);
    for (int ell = std::max(ell_min_, std::abs(m)); ell <= ell_max_; ++ell) {
        const teuk::Complex amplitude = coefficient(ell, m);
        for (std::size_t i = 0; i < z.size(); ++i)
            values[i] += amplitude * reduced_spin_weighted_spherical_harmonic(2, ell, m, z[i]);
    }
    return values;
}

std::vector<teuk::Complex> SchwarzschildPsi0Modes::exact_reduced_seed(
    int m, teuk::Real omega, const std::vector<teuk::Real>& z) const {
    std::vector<teuk::Complex> values(z.size(), teuk::zeroC);
    for (int ell = std::max(ell_min_, std::abs(m)); ell <= ell_max_; ++ell) {
        const teuk::Real starobinsky = teuk::Real(
            (ell + 2) * (ell + 1) * ell * (ell - 1));
        const teuk::Real denominator = starobinsky * starobinsky / 16 +
            9 * mass_ * mass_ * omega * omega;
        const teuk::Complex amplitude =
            teuk::Complex(0, 2) * omega * omega * omega * coefficient(ell, m) / denominator;
        for (std::size_t i = 0; i < z.size(); ++i)
            values[i] += amplitude * reduced_spin_weighted_spherical_harmonic(2, ell, m, z[i]);
    }
    return values;
}

} // namespace ghz::asymptotic
