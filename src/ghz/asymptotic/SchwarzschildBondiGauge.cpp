#include "ghz/asymptotic/SchwarzschildBondiGauge.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace ghz::asymptotic {
namespace {

int parity_sign(int value) { return (std::abs(value) % 2 == 0) ? 1 : -1; }

teuk::Complex phase(int m, teuk::Real omega, teuk::Real u, teuk::Real phi) {
    return std::exp(teuk::Complex(0, teuk::Real(m) * phi - omega * u));
}

} // namespace

HeldGaugeVectorMode SchwarzschildBondiGauge::held_mode(
    int ell, int m, teuk::Real omega) const {
    if (ell < 2 || std::abs(m) > ell)
        throw std::invalid_argument("Schwarzschild zeta mode requires ell >= 2 and |m| <= ell");

    const teuk::Real ell_r = ell;
    const teuk::Real starobinsky =
        (ell_r + 2) * (ell_r + 1) * ell_r * (ell_r - 1);
    const teuk::Real denominator = starobinsky * starobinsky / 16
        + 9 * psi0_.mass() * psi0_.mass() * omega * omega;
    const teuk::Complex psi = psi0_.coefficient(ell, m);
    const teuk::Complex electric_parity_psi = psi
        + teuk::Real(parity_sign(m)) * std::conj(psi0_.coefficient(ell, -m));
    const teuk::Real vector_factor =
        std::sqrt((ell_r + 2) * (ell_r - 1) / 2);
    const teuk::Real tensor_factor = std::sqrt(starobinsky / 4);

    HeldGaugeVectorMode result;
    result.m = teuk::Complex(0, -2) * omega * omega * omega
        * vector_factor * psi / denominator;
    result.l = omega * omega * tensor_factor * electric_parity_psi / denominator;
    result.n = -result.l * (ell_r * (ell_r + 1) / 2 + teuk::Real(0.5));
    return result;
}

HeldGaugeVectorGrid SchwarzschildBondiGauge::summed_reduced_held_mode(
    int m, teuk::Real omega, const std::vector<teuk::Real>& z) const {
    HeldGaugeVectorGrid result{
        std::vector<teuk::Complex>(z.size(), teuk::zeroC),
        std::vector<teuk::Complex>(z.size(), teuk::zeroC),
        std::vector<teuk::Complex>(z.size(), teuk::zeroC)};
    for (int ell = std::max({2, std::abs(m), psi0_.ell_min()});
         ell <= psi0_.ell_max(); ++ell) {
        const auto coefficient = held_mode(ell, m, omega);
        for (std::size_t i = 0; i < z.size(); ++i) {
            const teuk::Real scalar_harmonic =
                SchwarzschildPsi0Modes::reduced_spin_weighted_spherical_harmonic(
                    0, ell, m, z[i]);
            const teuk::Real vector_harmonic =
                SchwarzschildPsi0Modes::reduced_spin_weighted_spherical_harmonic(
                    1, ell, m, z[i]);
            result.l[i] += coefficient.l * scalar_harmonic;
            result.n[i] += coefficient.n * scalar_harmonic;
            result.m[i] += coefficient.m * vector_harmonic;
        }
    }
    return result;
}

NPGaugeVectorMode SchwarzschildBondiGauge::np_mode_at_radius(
    int ell, int m, teuk::Real omega, teuk::Real radius) const {
    if (radius <= teuk::Real(0))
        throw std::invalid_argument("gauge-vector evaluation radius must be positive");
    const HeldGaugeVectorMode held = held_mode(ell, m, omega);
    const HeldGaugeVectorMode reflected = held_mode(ell, -m, -omega);
    const teuk::Complex held_mbar = teuk::Real(parity_sign(m + 1))
        * std::conj(reflected.m);
    const teuk::Real eth_factor =
        std::sqrt(teuk::Real(ell * (ell + 1)) / 2);

    const teuk::Complex basis_l = held.n
        + teuk::Complex(0, radius * omega) * held.l
        - psi0_.mass() * held.l / radius;
    const teuk::Complex basis_n = held.l;
    const teuk::Complex basis_mbar = radius * held.m + eth_factor * held.l;
    const teuk::Complex basis_m = radius * held_mbar - eth_factor * held.l;
    return {basis_n, basis_l, -basis_mbar, -basis_m};
}

teuk::Complex SchwarzschildBondiGauge::u_dot_zeta_amplitude(
    int ell, int m, teuk::Real omega,
    const orbit::KerrCircularEquatorialOrbit& circular_orbit) const {
    const teuk::Real parameter_scale = std::max({
        teuk::Real(1), std::abs(psi0_.mass()), std::abs(circular_orbit.M()),
        std::abs(psi0_.orbital_radius()), std::abs(circular_orbit.radius())});
    const teuk::Real parameter_tolerance =
        teuk::Real(64) * std::numeric_limits<teuk::Real>::epsilon() * parameter_scale;
    if (std::abs(circular_orbit.a()) > parameter_tolerance
        || std::abs(circular_orbit.M() - psi0_.mass()) > parameter_tolerance
        || std::abs(circular_orbit.radius() - psi0_.orbital_radius())
            > parameter_tolerance) {
        throw std::invalid_argument(
            "Schwarzschild zeta data and circular trajectory parameters do not match");
    }
    const auto zeta = np_mode_at_radius(ell, m, omega, circular_orbit.radius());
    const auto velocity = circular_orbit.kinnersley_four_velocity();
    constexpr teuk::Real equator_z = 0;
    const teuk::Real y0 = SchwarzschildPsi0Modes::spin_weighted_spherical_harmonic(
        0, ell, m, equator_z);
    const teuk::Real y_plus = SchwarzschildPsi0Modes::spin_weighted_spherical_harmonic(
        1, ell, m, equator_z);
    const teuk::Real y_minus = SchwarzschildPsi0Modes::spin_weighted_spherical_harmonic(
        -1, ell, m, equator_z);
    return velocity.l * zeta.n * y0 + velocity.n * zeta.l * y0
        - velocity.m * zeta.mbar * y_minus
        - velocity.mbar * zeta.m * y_plus;
}

teuk::Complex SchwarzschildBondiGauge::lie_derivative_uu(
    int ell, int m, teuk::Real omega,
    const orbit::KerrCircularEquatorialOrbit& circular_orbit,
    teuk::Real outgoing_u, teuk::Real phi) const {
    const auto velocity = circular_orbit.outgoing_four_velocity();
    const teuk::Complex derivative_factor(
        0, velocity.u * (teuk::Real(m) * circular_orbit.Omega_phi() - omega));
    return teuk::Real(2) * derivative_factor
        * u_dot_zeta_amplitude(ell, m, omega, circular_orbit)
        * phase(m, omega, outgoing_u, phi);
}

} // namespace ghz::asymptotic
