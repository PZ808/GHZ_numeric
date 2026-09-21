#include "ghz/asymptotic/BondiHeldSeedSolve.hpp"

#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ghz::asymptotic {
namespace {

using Field = spectral::SpectralGHPVectorized;

Field make_field(std::size_t nz, int m, int p, int q, teuk::Real omega) {
    Field field(1, nz, {m, 0, 0}, GHPScalar<teuk::Complex>(teuk::zeroC, p, q), p, q);
    field.set_omega_mk(omega);
    return field;
}

teuk::Real infinity_norm(const Eigen::VectorXcd& value) {
    teuk::Real result = 0;
    for (Eigen::Index i = 0; i < value.size(); ++i)
        result = std::max(result, static_cast<teuk::Real>(std::abs(value(i))));
    return result;
}

} // namespace

BondiHeldSeedSolver::BondiHeldSeedSolver(
    const spectral::SpectralDiffer& differ,
    const KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>& held_operators,
    teuk::Real mass)
    : differ_(differ), held_operators_(held_operators), mass_(mass) {
    if (!(mass_ > 0)) throw std::invalid_argument("BondiHeldSeedSolver requires M > 0");
}

std::vector<teuk::Complex> BondiHeldSeedSolver::apply(
    const std::vector<teuk::Complex>& reduced_seed, int m, teuk::Real omega) const {
    const std::size_t nz = differ_.Nz();
    if (reduced_seed.size() != nz) throw std::invalid_argument("Bondi seed Nz mismatch");

    Field current = make_field(nz, m, 2, -2, omega);
    for (std::size_t i = 0; i < nz; ++i) current(0, i).value() = reduced_seed[i];

    int p = 2, q = -2;
    for (int count = 0; count < 4; ++count) {
        Field next = make_field(nz, m, p - 2, q, omega);
        held_operators_.edthBarHRed_inplace(current, next);
        current = std::move(next);
        p -= 2;
    }
    for (int count = 0; count < 4; ++count) {
        Field next = make_field(nz, m, p, q - 2, omega);
        held_operators_.edthHRed_inplace(current, next);
        current = std::move(next);
        q -= 2;
    }

    const teuk::Real mass_term = 9 * mass_ * mass_ * omega * omega;
    std::vector<teuk::Complex> result(nz);
    for (std::size_t i = 0; i < nz; ++i)
        result[i] = current(0, i).value() + mass_term * reduced_seed[i];
    return result;
}

BondiHeldSeedSolution BondiHeldSeedSolver::solve(
    const std::vector<teuk::Complex>& reduced_psi0, int m, teuk::Real omega) const {
    const Eigen::Index n = static_cast<Eigen::Index>(differ_.Nz());
    if (reduced_psi0.size() != static_cast<std::size_t>(n))
        throw std::invalid_argument("Bondi psi0 Nz mismatch");

    Eigen::MatrixXcd matrix(n, n);
    std::vector<teuk::Complex> basis(static_cast<std::size_t>(n), teuk::zeroC);
    for (Eigen::Index column = 0; column < n; ++column) {
        std::fill(basis.begin(), basis.end(), teuk::zeroC);
        basis[static_cast<std::size_t>(column)] = teuk::Complex(1, 0);
        const auto image = apply(basis, m, omega);
        for (Eigen::Index row = 0; row < n; ++row)
            matrix(row, column) = image[static_cast<std::size_t>(row)];
    }

    Eigen::VectorXcd rhs(n);
    const teuk::Complex rhs_factor(0, 2 * omega * omega * omega);
    for (Eigen::Index i = 0; i < n; ++i)
        rhs(i) = rhs_factor * reduced_psi0[static_cast<std::size_t>(i)];

    Eigen::VectorXd row_scale(n), column_scale(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        double maximum = 0;
        for (Eigen::Index j = 0; j < n; ++j) maximum = std::max(maximum, std::abs(matrix(i, j)));
        row_scale(i) = maximum > 0 ? 1 / maximum : 1;
    }
    Eigen::MatrixXcd equilibrated = row_scale.asDiagonal() * matrix;
    for (Eigen::Index j = 0; j < n; ++j) {
        double maximum = 0;
        for (Eigen::Index i = 0; i < n; ++i) maximum = std::max(maximum, std::abs(equilibrated(i, j)));
        column_scale(j) = maximum > 0 ? 1 / maximum : 1;
    }
    equilibrated *= column_scale.asDiagonal();
    const Eigen::VectorXcd equilibrated_rhs = row_scale.asDiagonal() * rhs;

    Eigen::FullPivLU<Eigen::MatrixXcd> lu(equilibrated);
    if (!lu.isInvertible()) throw std::runtime_error("Bondi Held collocation matrix is singular");
    Eigen::VectorXcd y = lu.solve(equilibrated_rhs);
    for (int refinement = 0; refinement < 3; ++refinement)
        y += lu.solve(equilibrated_rhs - equilibrated * y);
    const Eigen::VectorXcd solution = column_scale.asDiagonal() * y;

    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(equilibrated);
    const auto singular = svd.singularValues();
    const double condition = singular(singular.size() - 1) > 0
        ? singular(0) / singular(singular.size() - 1)
        : std::numeric_limits<double>::infinity();

    BondiHeldSeedSolution result;
    result.relative_residual = infinity_norm(matrix * solution - rhs) /
        std::max(infinity_norm(rhs), teuk::Real(1e-300));
    result.equilibrated_condition_number = condition;
    result.reduced_seed.resize(static_cast<std::size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i)
        result.reduced_seed[static_cast<std::size_t>(i)] = solution(i);
    return result;
}

} // namespace ghz::asymptotic
