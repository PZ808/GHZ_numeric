#include "ghz/asymptotic/BondiHeldMetricGrid.hpp"

#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"

#include <stdexcept>

namespace ghz::asymptotic {
namespace {

using Real = teuk::Real;
using Complex = teuk::Complex;
using Vector = std::vector<Complex>;
using Field = spectral::SpectralGHPVectorized;
constexpr Complex I(0, 1);

Vector zeros(std::size_t size) { return Vector(size, Complex(0)); }

Vector scale(const Vector& value, Complex factor) {
    Vector result(value.size());
    for (std::size_t i = 0; i < value.size(); ++i) result[i] = factor * value[i];
    return result;
}

Vector sum(const Vector& lhs, const Vector& rhs) {
    Vector result(lhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i) result[i] = lhs[i] + rhs[i];
    return result;
}

Vector difference(const Vector& lhs, const Vector& rhs) {
    return sum(lhs, scale(rhs, Complex(-1)));
}

MetricGridLaurent make_laurent(std::size_t nz) {
    return {{zeros(nz), zeros(nz), zeros(nz), zeros(nz)}};
}

MetricGridParts make_parts(std::size_t nz) {
    return {make_laurent(nz), make_laurent(nz), make_laurent(nz)};
}

void finish(MetricGridParts& parts) {
    for (std::size_t i = 0; i < 4; ++i)
        parts.adjusted.coefficient[i] = sum(parts.reconstructed.coefficient[i],
                                             parts.lie_zeta.coefficient[i]);
}

} // namespace

BondiHeldMetricReconstruction::BondiHeldMetricReconstruction(
    const spectral::SpectralDiffer& differ,
    const KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>& operators,
    Real mass)
    : differ_(differ), operators_(operators), mass_(mass) {
    if (!(mass_ > 0)) throw std::invalid_argument("Bondi metric requires M > 0");
}

Vector BondiHeldMetricReconstruction::up(
    const Vector& value, int spin, int m, Real omega) const {
    const std::size_t nz = differ_.Nz();
    if (value.size() != nz) throw std::invalid_argument("Held metric Nz mismatch");
    Field input(1, nz, {m, 0, 0}, GHPScalar<Complex>(teuk::zeroC, spin, -spin),
                spin, -spin);
    Field output(1, nz, {m, 0, 0}, GHPScalar<Complex>(teuk::zeroC, spin, -spin - 2),
                 spin, -spin - 2);
    input.set_omega_mk(omega);
    output.set_omega_mk(omega);
    for (std::size_t i = 0; i < nz; ++i) input(0, i).value() = value[i];
    operators_.edthHRed_inplace(input, output);
    Vector result(nz);
    for (std::size_t i = 0; i < nz; ++i) result[i] = output(0, i).value();
    return result;
}

Vector BondiHeldMetricReconstruction::down(
    const Vector& value, int spin, int m, Real omega) const {
    const std::size_t nz = differ_.Nz();
    if (value.size() != nz) throw std::invalid_argument("Held metric Nz mismatch");
    Field input(1, nz, {m, 0, 0}, GHPScalar<Complex>(teuk::zeroC, spin, -spin),
                spin, -spin);
    Field output(1, nz, {m, 0, 0}, GHPScalar<Complex>(teuk::zeroC, spin - 2, -spin),
                 spin - 2, -spin);
    input.set_omega_mk(omega);
    output.set_omega_mk(omega);
    for (std::size_t i = 0; i < nz; ++i) input(0, i).value() = value[i];
    operators_.edthBarHRed_inplace(input, output);
    Vector result(nz);
    for (std::size_t i = 0; i < nz; ++i) result[i] = output(0, i).value();
    return result;
}

BondiHeldMetricGrid BondiHeldMetricReconstruction::reconstruct(
    const Vector& f, const Vector& fbar, int m, Real omega) const {
    const std::size_t nz = differ_.Nz();
    if (f.size() != nz || fbar.size() != nz || omega == Real(0))
        throw std::invalid_argument("Held metric requires two Nz seeds and nonzero omega");
    const auto raise = [&](const Vector& v, int s) { return up(v, s, m, omega); };
    const auto lower = [&](const Vector& v, int s) { return down(v, s, m, omega); };
    const auto raise2 = [&](const Vector& v, int s) {
        return raise(raise(v, s), s + 1);
    };
    const auto lower2 = [&](const Vector& v, int s) {
        return lower(lower(v, s), s - 1);
    };
    const auto beta_half = [&](const Vector& v, int s) {
        return s == -2 ? scale(lower(raise(v, -2), -1), Complex(-1))
                       : scale(raise(lower(v, 2), 1), Complex(-1));
    };
    const auto lambda_quarter = [&](const Vector& v, int s) {
        return s == -2 ? lower2(raise2(v, -2), 0)
                       : raise2(lower2(v, 2), 0);
    };
    const auto potential = [&](const Vector& seed, int spin) {
        const Vector B = beta_half(seed, spin);
        const Vector Q = lambda_quarter(seed, spin);
        const Vector A = sum(beta_half(Q, spin), Q);
        std::array<Vector, 4> phi;
        phi[0] = scale(difference(A, scale(sum(B, seed),
                                              Real(3) * I * mass_ * omega)),
                       -Real(1) / (Real(6) * omega * omega));
        phi[1] = scale(difference(Q, scale(seed, Real(3) * I * mass_ * omega)),
                       Real(1) / (Real(2) * I * omega));
        phi[2] = B;
        phi[3] = scale(seed, I * omega);
        return phi;
    };
    const auto phi = potential(fbar, -2);
    const auto barphi = potential(f, 2);

    const Vector F0 = lower2(f, 2);
    const Vector Fbar0 = raise2(fbar, -2);
    const Vector zeta_l = scale(sum(F0, Fbar0), Real(1) / (Real(2) * I * omega));

    BondiHeldMetricGrid result{make_parts(nz), make_parts(nz), make_parts(nz)};
    const auto scalar_phi = [&](std::size_t i) {
        return sum(raise2(phi[i], -2), lower2(barphi[i], 2));
    };
    result.nn.reconstructed.coefficient[0] = scalar_phi(3);
    result.nn.reconstructed.coefficient[1] = scale(scalar_phi(2), Complex(-1));
    result.nn.reconstructed.coefficient[2] = scalar_phi(1);
    result.nn.reconstructed.coefficient[3] = scale(scalar_phi(0), Complex(-1));
    result.nm.reconstructed.coefficient[0] = scale(lower(barphi[3], 2), Complex(-1));
    result.nm.reconstructed.coefficient[2] = lower(barphi[1], 2);
    result.nm.reconstructed.coefficient[3] = scale(lower(barphi[0], 2), Complex(-2));
    result.mm.reconstructed.coefficient[1] = scale(barphi[2], Complex(2));
    result.mm.reconstructed.coefficient[2] = scale(barphi[1], Complex(-2));

    const Vector temp_f = sum(
        scale(lower(lower(lower(raise(f, 2), 3), 2), 1), Complex(-2)),
        scale(F0, Complex(4)));
    const Vector temp_fbar = sum(
        scale(raise(raise(raise(lower(fbar, -2), -3), -2), -1), Complex(-2)),
        scale(Fbar0, Complex(4)));
    result.nn.lie_zeta.coefficient[0] = scale(sum(F0, Fbar0), -I * omega);
    result.nn.lie_zeta.coefficient[1] = scale(sum(temp_f, temp_fbar), Complex(0.5));
    result.nn.lie_zeta.coefficient[2] = scale(sum(F0, Fbar0), Real(3) * mass_);
    result.nn.lie_zeta.coefficient[3] = scale(lower(raise(zeta_l, 0), 1), Real(2) * mass_);
    result.nm.lie_zeta.coefficient[0] = scale(lower(f, 2), I * omega);
    result.nm.lie_zeta.coefficient[2] = scale(lower(raise2(zeta_l, 0), 2), Complex(-1));
    result.nm.lie_zeta.coefficient[3] = scale(raise(zeta_l, 0), Real(2) * mass_);
    result.mm.lie_zeta.coefficient[1] = difference(
        scale(lower(raise(f, 2), 3), Complex(2)), scale(f, Complex(4)));
    result.mm.lie_zeta.coefficient[2] = scale(raise2(zeta_l, 0), Complex(2));
    finish(result.nn);
    finish(result.nm);
    finish(result.mm);
    return result;
}

} // namespace ghz::asymptotic
