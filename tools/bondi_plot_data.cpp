#include "ghz/asymptotic/BondiHeldSeedSolve.hpp"
#include "ghz/asymptotic/SchwarzschildBondiGauge.hpp"
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

using Real = teuk::Real;
using Complex = teuk::Complex;

Real relative_max_error(const std::vector<Complex>& actual,
                        const std::vector<Complex>& expected) {
    Real error = 0;
    Real scale = 0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        error = std::max(error, static_cast<Real>(std::abs(actual[i] - expected[i])));
        scale = std::max(scale, static_cast<Real>(std::abs(expected[i])));
    }
    return error / scale;
}

void write_redshift_data(const ghz::asymptotic::SchwarzschildPsi0Modes& psi0,
                         const std::filesystem::path& output) {
    KerrParams parameters(psi0.mass(), Real(0));
    KerrMetric metric(parameters);
    orbit::KerrCircularEquatorialOrbit trajectory(metric, psi0.orbital_radius(), +1);
    ghz::asymptotic::SchwarzschildBondiGauge gauge(psi0);
    const Real omega_orbit = trajectory.Omega_phi();
    constexpr Real detuning = Real(0.004);

    std::ofstream stream(output);
    if (!stream) throw std::runtime_error("cannot open redshift output: " + output.string());
    stream << std::setprecision(18)
           << "ell,m,omega_helical,omega_detuned,abs_u_dot_zeta_helical,"
              "lie_uu_helical_real,lie_uu_helical_imag,abs_lie_uu_helical,"
              "abs_lie_uu_detuned,abs_lie_uu_detuned_expected\n";
    Real max_helical = 0;
    Real max_detuned = 0;
    int count = 0;
    for (int ell = 2; ell <= std::min(8, psi0.ell_max()); ++ell) {
        for (int m = -ell; m <= ell; ++m) {
            if (m == 0) continue;
            const Real omega = Real(m) * omega_orbit;
            const Complex u_dot = gauge.u_dot_zeta_amplitude(ell, m, omega, trajectory);
            const Complex lie = gauge.lie_derivative_uu(
                ell, m, omega, trajectory, Real(0.31), Real(-0.27));
            const Complex lie_detuned = gauge.lie_derivative_uu(
                ell, m, omega + detuning, trajectory, Real(0.31), Real(-0.27));
            const Complex u_dot_detuned = gauge.u_dot_zeta_amplitude(
                ell, m, omega + detuning, trajectory);
            const Real detuned_expected = Real(2) * trajectory.gamma() * detuning
                * std::abs(u_dot_detuned);
            const Real detuned_error = std::abs(std::abs(lie_detuned) - detuned_expected);
            if (detuned_error > Real(2e-13) * std::max(Real(1), detuned_expected))
                throw std::runtime_error("detuned Lie derivative amplitude check failed");
            max_helical = std::max(max_helical, static_cast<Real>(std::abs(lie)));
            max_detuned = std::max(max_detuned, static_cast<Real>(std::abs(lie_detuned)));
            stream << ell << ',' << m << ',' << omega << ',' << omega + detuning << ','
                   << std::abs(u_dot) << ',' << lie.real() << ',' << lie.imag() << ','
                   << std::abs(lie) << ',' << std::abs(lie_detuned) << ','
                   << detuned_expected << '\n';
            ++count;
        }
    }
    std::cout << "redshift modes=" << count << " max helical |Lie_uu|=" << max_helical
              << " max detuned |Lie_uu|=" << max_detuned << '\n';
}

void write_m_mode_data(const ghz::asymptotic::SchwarzschildPsi0Modes& psi0,
                       const std::filesystem::path& grid_output,
                       const std::filesystem::path& summary_output) {
    constexpr std::size_t nz = 21;
    constexpr std::size_t nr = 4;
    spectral::SpectralDiffer differ(nz, nr);
    KerrParams parameters(psi0.mass(), Real(0));
    KerrMetric metric(parameters);
    CoordinateHelper coordinates(metric);
    KinnersleyTetrad<OutgoingCoords> tetrad(metric, coordinates);
    OutgoingCoords point(Real(0), psi0.orbital_radius(), Real(0), Real(0));
    auto held = ghp::build_held_fields_vectorized(tetrad, differ.lgl_nodes(), point).fields;
    const ghz::numeric::AffineMap1D radial_map({Real(2), Real(20)});
    KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer> operators(
        differ, held, tetrad, radial_map);
    ghz::asymptotic::BondiHeldSeedSolver solver(differ, operators, psi0.mass());

    std::ofstream grid(grid_output), summary(summary_output);
    if (!grid || !summary) throw std::runtime_error("cannot open m-mode output CSVs");
    grid << std::setprecision(18)
         << "m,z,psi0_reduced_real,psi0_reduced_imag,seed_numeric_real,"
            "seed_numeric_imag,seed_exact_real,seed_exact_imag,abs_seed_error\n";
    summary << std::setprecision(18)
            << "m,nz,omega,relative_residual,relative_seed_error,"
               "relative_exact_equation_residual,equilibrated_condition_number\n";
    const Real omega_orbit = std::sqrt(psi0.mass() /
        (psi0.orbital_radius() * psi0.orbital_radius() * psi0.orbital_radius()));
    for (int m : {2, 3, 4}) {
        const Real omega = Real(m) * omega_orbit;
        const auto input = psi0.summed_reduced_mode(m, differ.lgl_nodes());
        const auto exact = psi0.exact_reduced_seed(m, omega, differ.lgl_nodes());
        const auto solution = solver.solve(input, m, omega);
        const auto image = solver.apply(exact, m, omega);
        std::vector<Complex> rhs(nz);
        for (std::size_t i = 0; i < nz; ++i) {
            rhs[i] = Complex(0, Real(2) * omega * omega * omega) * input[i];
            grid << m << ',' << differ.lgl_nodes()[i] << ','
                 << input[i].real() << ',' << input[i].imag() << ','
                 << solution.reduced_seed[i].real() << ','
                 << solution.reduced_seed[i].imag() << ','
                 << exact[i].real() << ',' << exact[i].imag() << ','
                 << std::abs(solution.reduced_seed[i] - exact[i]) << '\n';
        }
        const Real seed_error = relative_max_error(solution.reduced_seed, exact);
        const Real exact_equation_error = relative_max_error(image, rhs);
        summary << m << ',' << nz << ',' << omega << ',' << solution.relative_residual
                << ',' << seed_error << ',' << exact_equation_error << ','
                << solution.equilibrated_condition_number << '\n';
        std::cout << "m=" << m << " collocation residual=" << solution.relative_residual
                  << " seed error=" << seed_error
                  << " exact equation residual=" << exact_equation_error << '\n';
        if (!(solution.relative_residual < Real(2e-6)
              && seed_error < Real(2e-7)
              && exact_equation_error < Real(5e-6)))
            throw std::runtime_error("m-mode diagnostics failed an existing solve tolerance");
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("usage: bondi_plot_data psi0.csv output_directory");
        const auto psi0 = ghz::asymptotic::SchwarzschildPsi0Modes::load_csv(argv[1]);
        const std::filesystem::path output(argv[2]);
        std::filesystem::create_directories(output);
        std::cout << std::setprecision(17);
        write_redshift_data(psi0, output / "redshift_modes.csv");
        write_m_mode_data(psi0, output / "m_mode_grid.csv", output / "m_mode_summary.csv");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
