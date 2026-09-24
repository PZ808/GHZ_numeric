#include "ghz/asymptotic/BondiHeldSeedSolve.hpp"
#include "ghz/asymptotic/BondiHeldMetricGrid.hpp"
#include "ghz/asymptotic/SchwarzschildBondiMetric.hpp"
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
#include <fstream>
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

teuk::Complex at_radius(const ghz::asymptotic::MetricGridLaurent& series,
                        std::size_t node, teuk::Real radius) {
    return series.coefficient[0][node] * radius
        + series.coefficient[1][node]
        + series.coefficient[2][node] / radius
        + series.coefficient[3][node] / (radius * radius);
}

void write_radial_rows(std::ostream& stream, int m, teuk::Real z,
                       teuk::Real pole, const char* name,
                       const ghz::asymptotic::MetricGridParts& parts,
                       std::size_t node) {
    for (int i = 0; i <= 60; ++i) {
        const teuk::Real radius = teuk::Real(20)
            * std::pow(teuk::Real(100), teuk::Real(i) / teuk::Real(60));
        const auto reconstructed = pole * at_radius(parts.reconstructed, node, radius);
        const auto lie = pole * at_radius(parts.lie_zeta, node, radius);
        const auto adjusted = reconstructed + lie;
        stream << m << ',' << name << ',' << z << ',' << radius << ','
               << std::abs(reconstructed) << ',' << std::abs(lie) << ','
               << std::abs(adjusted) << ','
               << std::abs(pole * parts.adjusted.coefficient[2][node]) << '\n';
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3)
            throw std::runtime_error("usage: bondi_schwarzschild_data_collocation_tests psi0.csv [spectral_plot.csv]");
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
        const auto negative_input = modes.summed_reduced_mode(-m, differ.lgl_nodes());
        const auto negative_exact = modes.exact_reduced_seed(-m, -omega, differ.lgl_nodes());
        const auto negative_solution = solver.solve(negative_input, -m, -omega);
        const teuk::Real negative_error = relative_max_error(
            negative_solution.reduced_seed, negative_exact);
        std::cout << "m=" << -m << " negative-mode seed error = " << negative_error << '\n';
        if (!(negative_solution.relative_residual < teuk::Real(2e-6)
              && negative_error < teuk::Real(2e-7)))
            throw std::runtime_error("negative-m Held collocation disagrees with Schwarzschild inversion");

        ghz::asymptotic::BondiHeldMetricReconstruction reconstruction(differ, operators, mass);
        ghz::asymptotic::SchwarzschildBondiMetric modal_reference(modes);
        std::ofstream plot_stream;
        if (argc == 3) {
            plot_stream.open(argv[2]);
            if (!plot_stream) throw std::runtime_error("cannot open spectral plot CSV");
            plot_stream << std::setprecision(18)
                        << "m,component,z,r,abs_reconstructed,abs_lie_zeta,abs_direct_sum,"
                           "abs_adjusted_coefficient_1_over_r\n";
        }
        for (int metric_m : {2, 3}) {
            const teuk::Real metric_omega = teuk::Real(metric_m) / std::sqrt(teuk::Real(1000));
            const auto positive_solve = solver.solve(
                modes.summed_reduced_mode(metric_m, differ.lgl_nodes()),
                metric_m, metric_omega);
            const auto negative_solve = solver.solve(
                modes.summed_reduced_mode(-metric_m, differ.lgl_nodes()),
                -metric_m, -metric_omega);
            const auto& positive_seed = positive_solve.reduced_seed;
            const auto& reflected_seed = negative_solve.reduced_seed;
            std::cout << "m=" << metric_m << " seed errors +,- = "
                      << relative_max_error(positive_seed,
                          modes.exact_reduced_seed(metric_m, metric_omega, differ.lgl_nodes()))
                      << ", " << relative_max_error(reflected_seed,
                          modes.exact_reduced_seed(-metric_m, -metric_omega, differ.lgl_nodes()))
                      << '\n';
            std::vector<teuk::Complex> reduced_fbar(nz);
            // The reduced spin-harmonic reflection already carries (-1)^m;
            // conjugating the solved -m grid is the complete +m partner.
            for (std::size_t i = 0; i < nz; ++i)
                reduced_fbar[i] = std::conj(reflected_seed[i]);
            const auto grid = reconstruction.reconstruct(
                positive_seed, reduced_fbar, metric_m, metric_omega);
            teuk::Real largest_reference_difference = 0;
            teuk::Real largest_cancellation = 0;
            std::string largest_label;
            for (std::size_t i = 1; i + 1 < nz; ++i) {
                const teuk::Real z = differ.lgl_nodes()[i];
                const auto reference = modal_reference.summed_at_z(metric_m, metric_omega, z);
                const struct {
                    const ghz::asymptotic::MetricGridParts* grid;
                    const ghz::asymptotic::MetricComponentParts* reference;
                    int spin;
                } components[]{
                    {&grid.nn, &reference.nn, 0},
                    {&grid.nm, &reference.nm, 1},
                    {&grid.mm, &reference.mm, 2}};
                for (const auto& component : components) {
                    const teuk::Real pole = modes.pole_factor(metric_m, component.spin, z);
                    for (std::size_t power = 0; power < 4; ++power) {
                        const auto calculated = pole * component.grid->reconstructed.coefficient[power][i];
                        const auto expected = component.reference->reconstructed.coefficient[power];
                        const teuk::Real difference = std::abs(calculated - expected)
                            / std::max(teuk::Real(1e-8), std::abs(expected));
                        if (difference > largest_reference_difference) {
                            largest_reference_difference = difference;
                            largest_label = std::to_string(component.spin) + ":"
                                + std::to_string(power) + ":" + std::to_string(i);
                        }
                        if (power < 2) {
                            const teuk::Real scale = std::max(
                                std::abs(component.grid->reconstructed.coefficient[power][i]),
                                std::abs(component.grid->lie_zeta.coefficient[power][i]));
                            if (scale > teuk::Real(1e-13))
                                largest_cancellation = std::max(
                                    largest_cancellation,
                                    static_cast<teuk::Real>(
                                        std::abs(component.grid->adjusted.coefficient[power][i])
                                        / scale));
                        }
                    }
                }
            }
            std::cout << "m=" << metric_m
                      << " spectral metric/reference difference=" << largest_reference_difference
                      << " at " << largest_label
                      << " growing/constant cancellation=" << largest_cancellation << '\n';
            if (!(largest_reference_difference < teuk::Real(1e-5)
                  && largest_cancellation < teuk::Real(1e-8)))
                throw std::runtime_error("spectral Held metric disagrees with modal benchmark");

            const auto nearest_node = std::min_element(
                differ.lgl_nodes().begin() + 1, differ.lgl_nodes().end() - 1,
                [](teuk::Real a, teuk::Real b) {
                    return std::abs(a - teuk::Real(0.37))
                        < std::abs(b - teuk::Real(0.37));
                });
            const std::size_t node = std::distance(differ.lgl_nodes().begin(), nearest_node);
            const teuk::Real z = *nearest_node;
            const std::pair<const char*, const ghz::asymptotic::MetricGridParts*> grid_components[]{
                {"nn", &grid.nn}, {"nm", &grid.nm}, {"mm", &grid.mm}};
            for (std::size_t spin = 0; spin < 3; ++spin) {
                const auto& [name, parts] = grid_components[spin];
                const teuk::Real pole = modes.pole_factor(metric_m, static_cast<int>(spin), z);
                const auto value = [&](teuk::Real radius) {
                    return pole * (at_radius(parts->reconstructed, node, radius)
                                   + at_radius(parts->lie_zeta, node, radius));
                };
                const teuk::Real slope = std::log(std::abs(value(teuk::Real(2000)))
                    / std::abs(value(teuk::Real(1000)))) / std::log(teuk::Real(2));
                std::cout << "m=" << metric_m << " h_" << name
                          << " spectral direct-sum slope=" << slope
                          << " at z=" << z << '\n';
                if (!(std::abs(slope + teuk::Real(1)) < teuk::Real(0.01)))
                    throw std::runtime_error("spectral Held metric does not fall as 1/r");
                if (plot_stream)
                    write_radial_rows(plot_stream, metric_m, z, pole, name, *parts, node);
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
