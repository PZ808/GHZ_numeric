#include "ghz/asymptotic/BondiHeldSeedSolve.hpp"
#include "ghz/asymptotic/SchwarzschildPsi0Modes.hpp"
#include "ghz/geom/Coords.hpp"
#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/KerrParams.hpp"
#include "ghz/geom/KinnersleyTetrad.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/spectral/SpectralCoordinateMaps.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

teuk::Real relative_max_error(const std::vector<teuk::Complex>& actual,
                              const std::vector<teuk::Complex>& expected,
                              const std::vector<teuk::Real>& weights = {}) {
    teuk::Real error = 0;
    teuk::Real scale = 0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        const teuk::Real weight = weights.empty() ? teuk::Real(1) : weights[i];
        error = std::max(error, static_cast<teuk::Real>(std::abs(weight * (actual[i] - expected[i]))));
        scale = std::max(scale, static_cast<teuk::Real>(std::abs(weight * expected[i])));
    }
    return error / scale;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("expected path to the psi0 coefficient CSV");
        constexpr std::size_t nz = 21;
        constexpr std::size_t nr = 4;
        constexpr int m = 2;
        constexpr teuk::Real mass = 1;
        constexpr teuk::Real spin = 0;
        const teuk::Real omega = teuk::Real(m) / std::sqrt(teuk::Real(1000));

        const auto modes = ghz::asymptotic::SchwarzschildPsi0Modes::load_csv(argv[1]);
        spectral::SpectralDiffer differ(nz, nr);
        KerrParams params(mass, spin);
        KerrMetric metric(params);
        CoordinateHelper coordinates(metric);
        KinnersleyTetrad<OutgoingCoords> tetrad(metric, coordinates);
        OutgoingCoords point(teuk::Real(0), teuk::Real(10), teuk::Real(0), teuk::Real(0));
        auto held = ghp::build_held_fields_vectorized(tetrad, differ.lgl_nodes(), point).fields;
        const ghz::numeric::Domain domain{teuk::Real(2), teuk::Real(20)};
        const ghz::numeric::AffineMap1D radial_map(domain);
        KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer> operators(
            differ, held, tetrad, radial_map);
        ghz::asymptotic::BondiHeldSeedSolver solver(differ, operators, mass);

        const auto psi0 = modes.summed_reduced_mode(m, differ.lgl_nodes());
        const auto exact = modes.exact_reduced_seed(m, omega, differ.lgl_nodes());
        const auto result = solver.solve(psi0, m, omega);
        const auto exact_image = solver.apply(exact, m, omega);

        std::vector<teuk::Complex> rhs(nz);
        std::vector<teuk::Real> pole(nz);
        for (std::size_t i = 0; i < nz; ++i) {
            rhs[i] = teuk::Complex(0, 2 * omega * omega * omega) * psi0[i];
            pole[i] = modes.pole_factor(m, 2, differ.lgl_nodes()[i]);
        }

        const teuk::Real solution_error = relative_max_error(result.reduced_seed, exact, pole);
        const teuk::Real exact_equation_error = relative_max_error(exact_image, rhs);
        std::cout << std::setprecision(17)
                  << "m=" << m << " Nz=" << nz << " omega=" << omega << '\n'
                  << "equilibrated condition number = "
                  << result.equilibrated_condition_number << '\n'
                  << "collocation residual = " << result.relative_residual << '\n'
                  << "raw seed error = " << solution_error << '\n'
                  << "exact lm-sum equation residual = " << exact_equation_error << '\n';

        if (!(result.relative_residual < teuk::Real(2e-6)))
            throw std::runtime_error("collocation residual is too large");
        if (!(solution_error < teuk::Real(2e-7)))
            throw std::runtime_error("collocation seed disagrees with the exact Schwarzschild inversion");
        if (!(exact_equation_error < teuk::Real(5e-6)))
            throw std::runtime_error("Held operator chain disagrees with the exact lm sum");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
