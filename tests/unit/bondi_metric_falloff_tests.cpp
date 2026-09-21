#include "ghz/asymptotic/SchwarzschildBondiMetric.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

using Real = teuk::Real;
using ghz::asymptotic::MetricComponentParts;

Real cancellation_error(const MetricComponentParts& value) {
    Real worst = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        const Real scale = std::max(std::abs(value.reconstructed.coefficient[i]),
                                    std::abs(value.lie_zeta.coefficient[i]));
        if (scale > 0)
            worst = std::max(worst, static_cast<Real>(
                std::abs(value.adjusted.coefficient[i]) / scale));
    }
    return worst;
}

void write_plot_data(std::ostream& stream,
                     const ghz::asymptotic::SchwarzschildBondiMetric& metric,
                     Real omega_orbit) {
    stream << std::setprecision(18)
           << "m,component,z,r,abs_reconstructed,abs_lie_zeta,abs_direct_sum,abs_adjusted,"
              "abs_adjusted_coefficient_1_over_r\n";
    constexpr Real z = Real(0.37);
    for (int m : {2, 3}) {
        const auto summed = metric.summed_at_z(m, Real(m) * omega_orbit, z);
        const std::pair<const char*, const MetricComponentParts*> components[]{
            {"nn", &summed.nn}, {"nm", &summed.nm}, {"mm", &summed.mm}};
        for (const auto& [name, value] : components) {
            for (int i = 0; i <= 60; ++i) {
                const Real radius = Real(20) * std::pow(Real(100), Real(i) / Real(60));
                stream << m << ',' << name << ',' << z << ',' << radius << ','
                       << std::abs(value->reconstructed.at(radius)) << ','
                       << std::abs(value->lie_zeta.at(radius)) << ','
                       << std::abs(value->reconstructed.at(radius)
                                   + value->lie_zeta.at(radius)) << ','
                       << std::abs(value->adjusted.at(radius)) << ','
                       << std::abs(value->adjusted.coefficient[2]) << '\n';
            }
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3)
            throw std::invalid_argument("usage: bondi_metric_falloff_tests psi0.csv [plot.csv]");
        const auto psi0 = ghz::asymptotic::SchwarzschildPsi0Modes::load_csv(argv[1]);
        ghz::asymptotic::SchwarzschildBondiMetric metric(psi0);
        const Real omega_orbit = std::sqrt(psi0.mass() /
            (psi0.orbital_radius() * psi0.orbital_radius() * psi0.orbital_radius()));
        Real max_cancellation_error = 0;
        Real max_surviving_inverse_radius = 0;
        int checked = 0;
        for (int ell = psi0.ell_min(); ell <= psi0.ell_max(); ++ell) {
            for (int m = -ell; m <= ell; ++m) {
                if (m == 0) continue;
                const auto result = metric.mode(ell, m, Real(m) * omega_orbit);
                for (const auto* value : {&result.nn, &result.nm, &result.mm}) {
                    max_cancellation_error = std::max(
                        max_cancellation_error, cancellation_error(*value));
                    max_surviving_inverse_radius = std::max(
                        max_surviving_inverse_radius,
                        static_cast<Real>(std::abs(value->adjusted.coefficient[2])));
                    ++checked;
                }
            }
        }
        std::cout << std::setprecision(17)
                  << "checked radiative (ell,m) components=" << checked << '\n'
                  << "max relative r and constant cancellation="
                  << max_cancellation_error << '\n'
                  << "max surviving |1/r coefficient|="
                  << max_surviving_inverse_radius << '\n';
        if (!(max_cancellation_error < Real(5e-14)))
            throw std::runtime_error("growing or constant metric term failed to cancel");
        if (!(max_surviving_inverse_radius > Real(1e-8)))
            throw std::runtime_error("falloff check has no nonzero 1/r signal");
        for (int m : {2, 3}) {
            const auto summed = metric.summed_at_z(m, Real(m) * omega_orbit, Real(0.37));
            for (const auto* value : {&summed.nn, &summed.nm, &summed.mm}) {
                const Real first = std::abs(value->reconstructed.at(Real(1000))
                    + value->lie_zeta.at(Real(1000)));
                const Real second = std::abs(value->reconstructed.at(Real(2000))
                    + value->lie_zeta.at(Real(2000)));
                const Real slope = std::log(second / first) / std::log(Real(2));
                if (!(std::abs(slope + Real(1)) < Real(0.01)))
                    throw std::runtime_error("direct finite-radius metric does not fall as 1/r");
            }
        }
        if (argc == 3) {
            std::ofstream stream(argv[2]);
            if (!stream) throw std::runtime_error("cannot open metric plot CSV");
            write_plot_data(stream, metric, omega_orbit);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
