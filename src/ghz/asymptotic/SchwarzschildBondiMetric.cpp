#include "ghz/asymptotic/SchwarzschildBondiMetric.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ghz::asymptotic {
namespace {

using Real = teuk::Real;
using Complex = teuk::Complex;
constexpr Complex I(0, 1);

int parity(int m) { return std::abs(m) % 2 == 0 ? 1 : -1; }

struct Potential {
    std::array<Complex, 4> phi; // Phi^0, Phi^1, Phi^2, Phi^3
    Complex f;
    Complex zeta_l;
};

Potential potential(const SchwarzschildPsi0Modes& psi0, int ell, int m,
                    Real omega) {
    const Real L = Real(ell) * Real(ell + 1);
    const Real beta = Real(ell + 2) * Real(ell - 1);
    const Real Lambda = L * beta;
    const Real D = Lambda * Lambda / Real(16)
        + Real(9) * psi0.mass() * psi0.mass() * omega * omega;
    const Complex psi = psi0.coefficient(ell, m);
    const Complex reflected_psi = Real(parity(m))
        * std::conj(psi0.coefficient(ell, -m));
    const Real omega2 = omega * omega;
    const Real omega3 = omega2 * omega;
    const Real omega4 = omega2 * omega2;
    Potential result;
    result.f = Real(2) * I * omega3 * psi / D;
    result.phi[0] = (-I * omega / (Real(3) * D))
        * (Lambda * L / Real(8) - Real(3) * I * psi0.mass() * omega * L / Real(2))
        * reflected_psi;
    result.phi[1] = (omega2 / D)
        * (Lambda / Real(4) - Real(3) * I * psi0.mass() * omega)
        * reflected_psi;
    result.phi[2] = I * omega3 * beta * reflected_psi / D;
    result.phi[3] = -Real(2) * omega4 * reflected_psi / D;
    result.zeta_l = omega2 * std::sqrt(Lambda / Real(4))
        * (psi + reflected_psi) / D;
    return result;
}

RadialLaurent add(const RadialLaurent& lhs, const RadialLaurent& rhs) {
    RadialLaurent result;
    for (std::size_t i = 0; i < result.coefficient.size(); ++i)
        result.coefficient[i] = lhs.coefficient[i] + rhs.coefficient[i];
    return result;
}

MetricComponentParts parts(RadialLaurent reconstructed, RadialLaurent lie) {
    return {reconstructed, lie, add(reconstructed, lie)};
}

} // namespace

Complex RadialLaurent::at(Real radius) const {
    if (!(radius > 0)) throw std::invalid_argument("metric evaluation requires r > 0");
    return coefficient[0] * radius + coefficient[1]
        + coefficient[2] / radius + coefficient[3] / (radius * radius);
}

SchwarzschildBondiMetricMode SchwarzschildBondiMetric::mode(
    int ell, int m, Real omega) const {
    if (ell < 2 || std::abs(m) > ell)
        throw std::invalid_argument("Bondi metric mode requires ell >= 2 and |m| <= ell");
    const Potential direct = potential(psi0_, ell, m, omega);
    const Potential reflected = potential(psi0_, ell, -m, -omega);
    std::array<Complex, 4> bar_phi;
    for (std::size_t i = 0; i < bar_phi.size(); ++i)
        bar_phi[i] = Real(parity(m)) * std::conj(reflected.phi[i]);

    const Real L = Real(ell) * Real(ell + 1);
    const Real beta = Real(ell + 2) * Real(ell - 1);
    const Real c = std::sqrt(L * beta) / Real(2); // two Held eth operators
    const Real up0 = std::sqrt(L / Real(2));
    const Real up1 = std::sqrt(beta / Real(2));
    const Real mass = psi0_.mass();

    // 2 Re S-dagger(Phi^Adj), with Phi^Adj = Phi0-r Phi1+r^2 Phi2-r^3 Phi3.
    RadialLaurent s_nn{{
        c * (direct.phi[3] + bar_phi[3]),
        -c * (direct.phi[2] + bar_phi[2]),
        c * (direct.phi[1] + bar_phi[1]),
        -c * (direct.phi[0] + bar_phi[0])}};
    RadialLaurent s_nm{{
        up1 * bar_phi[3], Complex(0), -up1 * bar_phi[1],
        Real(2) * up1 * bar_phi[0]}};
    RadialLaurent s_mm{{
        Complex(0), Real(2) * bar_phi[2], -Real(2) * bar_phi[1], Complex(0)}};

    // Schwarzschild reduction of LieZetagComps in BondiGauge_Adjusted.nb.
    // The nn expression is a real spin-zero modal projection of Temp.
    const auto nn_temp = [&](const Potential& value, Real mode_omega) {
        return RadialLaurent{{
            -Real(2) * I * mode_omega * c * value.f,
            beta * c * value.f,
            Real(6) * mass * c * value.f,
            -mass * L * value.zeta_l}};
    };
    const auto temp = nn_temp(direct, omega);
    const auto reflected_temp = nn_temp(reflected, -omega);
    RadialLaurent lie_nn;
    for (std::size_t i = 0; i < lie_nn.coefficient.size(); ++i)
        lie_nn.coefficient[i] = (temp.coefficient[i] + Real(parity(m))
            * std::conj(reflected_temp.coefficient[i])) / Real(2);
    RadialLaurent lie_nm{{
        -I * omega * up1 * direct.f,
        Complex(0),
        beta * up0 * direct.zeta_l / Real(2),
        Real(2) * mass * up0 * direct.zeta_l}};
    RadialLaurent lie_mm{{
        Complex(0), -beta * direct.f,
        Real(2) * c * direct.zeta_l, Complex(0)}};

    return {parts(s_nn, lie_nn), parts(s_nm, lie_nm), parts(s_mm, lie_mm)};
}

SchwarzschildBondiMetricMode SchwarzschildBondiMetric::summed_at_z(
    int m, Real omega, Real z) const {
    SchwarzschildBondiMetricMode result;
    const int first_ell = std::max({2, std::abs(m), psi0_.ell_min()});
    for (int ell = first_ell; ell <= psi0_.ell_max(); ++ell) {
        const auto current = mode(ell, m, omega);
        const std::array<Real, 3> harmonics{
            SchwarzschildPsi0Modes::spin_weighted_spherical_harmonic(0, ell, m, z),
            SchwarzschildPsi0Modes::spin_weighted_spherical_harmonic(1, ell, m, z),
            SchwarzschildPsi0Modes::spin_weighted_spherical_harmonic(2, ell, m, z)};
        auto accumulate = [&](MetricComponentParts& target,
                              const MetricComponentParts& source, Real harmonic) {
            for (std::size_t i = 0; i < 4; ++i) {
                target.reconstructed.coefficient[i] += harmonic * source.reconstructed.coefficient[i];
                target.lie_zeta.coefficient[i] += harmonic * source.lie_zeta.coefficient[i];
                target.adjusted.coefficient[i] += harmonic * source.adjusted.coefficient[i];
            }
        };
        accumulate(result.nn, current.nn, harmonics[0]);
        accumulate(result.nm, current.nm, harmonics[1]);
        accumulate(result.mm, current.mm, harmonics[2]);
    }
    return result;
}

} // namespace ghz::asymptotic
