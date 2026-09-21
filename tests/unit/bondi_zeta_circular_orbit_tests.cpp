#include "ghz/asymptotic/SchwarzschildBondiGauge.hpp"
#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/KerrParams.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using Real = teuk::Real;
using Complex = teuk::Complex;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool close(Real actual, Real expected, Real tolerance) {
    return std::abs(actual - expected) <= tolerance;
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("expected path to the psi0 coefficient CSV");
        constexpr Real mass = 1;
        constexpr Real radius = 10;
        constexpr Real tolerance = 2e-15L;

        KerrParams parameters(mass, Real(0));
        KerrMetric metric(parameters);
        orbit::KerrCircularEquatorialOrbit trajectory(metric, radius, +1);

        const Real expected_omega = std::sqrt(mass / (radius * radius * radius));
        const Real expected_gamma = std::sqrt(radius / (radius - Real(3) * mass));
        require(close(trajectory.Omega_phi(), expected_omega, tolerance),
                "circular Schwarzschild orbital frequency is incorrect");
        require(close(trajectory.gamma(), expected_gamma, tolerance),
                "circular Schwarzschild gamma is incorrect");
        const auto frequencies = trajectory.source_frequencies();
        require(close(frequencies.Omega_phi, expected_omega, tolerance)
                    && frequencies.Omega_r == Real(0) && frequencies.Omega_z == Real(0),
                "circular source frequencies are assigned to the wrong fields");

        const auto velocity = trajectory.kinnersley_four_velocity();
        const Complex velocity_norm = Real(2) * velocity.l * velocity.n
            - Real(2) * velocity.m * velocity.mbar;
        require(std::abs(velocity_norm - Complex(1, 0)) < tolerance,
                "mostly-minus Kinnersley velocity does not have unit norm");

        const Real inverse_sqrt_four_pi = Real(1) /
            std::sqrt(Real(4) * std::acos(Real(-1)));
        require(close(
                    ghz::asymptotic::SchwarzschildPsi0Modes::
                        spin_weighted_spherical_harmonic(0, 0, 0, Real(0.37)),
                    inverse_sqrt_four_pi, tolerance),
                "scalar Y00 convention is incorrect");
        require(close(
                    ghz::asymptotic::SchwarzschildPsi0Modes::
                        spin_weighted_spherical_harmonic(0, 1, 0, Real(0.37)),
                    std::sqrt(Real(3)) * Real(0.37) * inverse_sqrt_four_pi,
                    tolerance),
                "scalar Y10 convention is incorrect");
        const Real spin_one_equator =
            -std::sqrt(Real(3) / (Real(16) * std::acos(Real(-1))));
        require(close(
                    ghz::asymptotic::SchwarzschildPsi0Modes::
                        spin_weighted_spherical_harmonic(1, 1, 1, Real(0)),
                    spin_one_equator, tolerance)
                    && close(
                        ghz::asymptotic::SchwarzschildPsi0Modes::
                            spin_weighted_spherical_harmonic(-1, 1, 1, Real(0)),
                        spin_one_equator, tolerance),
                "spin-one spherical-harmonic convention is incorrect");

        // Cross-check the general Wigner-d implementation against the older
        // Jacobi implementation used by the psi0 collocation solve.
        for (int ell = 2; ell <= 20; ++ell) {
            for (int m = 2; m <= ell; ++m) {
                for (const Real z : {Real(-1), Real(-0.37), Real(0), Real(0.44), Real(1)}) {
                    const Real general =
                        ghz::asymptotic::SchwarzschildPsi0Modes::
                            reduced_spin_weighted_spherical_harmonic(2, ell, m, z);
                    const Real jacobi =
                        ghz::asymptotic::SchwarzschildPsi0Modes::
                            reduced_spin2_spherical_harmonic(ell, m, z);
                    const Real harmonic_scale = std::max(Real(1), std::abs(jacobi));
                    if (std::abs(general - jacobi) >= Real(2e-13L) * harmonic_scale) {
                        std::cerr << "harmonic mismatch ell=" << ell << " m=" << m
                                  << " z=" << z << " general=" << general
                                  << " jacobi=" << jacobi << '\n';
                        throw std::runtime_error(
                            "general spin-weighted harmonic convention disagrees with psi0");
                    }
                }
            }
        }

        const auto psi0 = ghz::asymptotic::SchwarzschildPsi0Modes::load_csv(argv[1]);
        ghz::asymptotic::SchwarzschildBondiGauge gauge(psi0);
        const std::vector<std::pair<int, int>> notebook_modes{
            {2, -2}, {2, -1}, {2, 1}, {2, 2},
            {3, -3}, {3, -2}, {3, -1}, {3, 1}, {3, 2}, {3, 3}};

        Real max_u_dot_zeta = 0;
        Real max_lie_uu = 0;
        for (const auto [ell, m] : notebook_modes) {
            const Real omega = Real(m) * frequencies.Omega_phi;
            const Complex u_dot_zeta =
                gauge.u_dot_zeta_amplitude(ell, m, omega, trajectory);
            const Complex lie_uu = gauge.lie_derivative_uu(
                ell, m, omega, trajectory, Real(0.31), Real(-0.27));
            max_u_dot_zeta = std::max(max_u_dot_zeta, static_cast<Real>(std::abs(u_dot_zeta)));
            max_lie_uu = std::max(max_lie_uu, static_cast<Real>(std::abs(lie_uu)));
            require(std::abs(lie_uu) < Real(2e-28L),
                    "(L_zeta g)_uu is not zero for a helical circular mode");
        }

        // Move one mode away from omega=m Omega and compare the analytic
        // modal Lie derivative with a centered derivative along the orbit.
        constexpr int off_ell = 2;
        constexpr int off_m = 2;
        const Real off_omega = Real(off_m) * frequencies.Omega_phi + Real(0.004);
        const Complex off_amplitude =
            gauge.u_dot_zeta_amplitude(off_ell, off_m, off_omega, trajectory);
        const Real step = Real(1e-4L);
        const auto pullback = [&](Real proper_time) {
            const Real u = trajectory.gamma() * proper_time;
            const Real phi = trajectory.gamma() * trajectory.Omega_phi() * proper_time;
            return off_amplitude * std::exp(
                Complex(0, Real(off_m) * phi - off_omega * u));
        };
        const Complex finite_difference =
            (pullback(step) - pullback(-step)) / step;
        const Complex analytic_off_helical = gauge.lie_derivative_uu(
            off_ell, off_m, off_omega, trajectory);
        require(std::abs(analytic_off_helical) > Real(1e-8L),
                "off-helical Lie-derivative check did not exercise a nonzero mode");
        require(std::abs(analytic_off_helical - finite_difference)
                    < Real(2e-10L) * std::abs(analytic_off_helical),
                "modal Lie derivative disagrees with the worldline finite difference");

        // Exercise the ell-summed, pole-factorized z-grid interface used by
        // the spectral Held operators, including both endpoints.
        const std::vector<Real> z_grid{Real(-1), Real(-0.5), Real(0), Real(0.5), Real(1)};
        const auto reduced = gauge.summed_reduced_held_mode(
            2, Real(2) * frequencies.Omega_phi, z_grid);
        require(reduced.l.size() == z_grid.size()
                    && reduced.n.size() == z_grid.size()
                    && reduced.m.size() == z_grid.size(),
                "summed zeta grid has the wrong extent");
        for (std::size_t i = 0; i < z_grid.size(); ++i) {
            require(std::isfinite(reduced.l[i].real()) && std::isfinite(reduced.l[i].imag())
                        && std::isfinite(reduced.n[i].real()) && std::isfinite(reduced.n[i].imag())
                        && std::isfinite(reduced.m[i].real()) && std::isfinite(reduced.m[i].imag()),
                    "summed reduced zeta mode is singular on the z grid");
        }

        std::cout.precision(17);
        std::cout << "Omega_phi=" << trajectory.Omega_phi()
                  << " gamma=" << trajectory.gamma() << '\n'
                  << "max |u.zeta| over notebook modes = " << max_u_dot_zeta << '\n'
                  << "max |(L_zeta g)_uu| over notebook modes = " << max_lie_uu << '\n'
                  << "off-helical |(L_zeta g)_uu| = "
                  << std::abs(analytic_off_helical) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
