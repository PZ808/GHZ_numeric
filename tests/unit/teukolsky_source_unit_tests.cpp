//
// teukolsky_source_unit_tests.cpp
//

#include "ghz/source/AngularProjection.hpp"
#include "ghz/source/TeukolskySource.hpp"

#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/KinnersleyTetrad.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

    using ghz::source::AngularModeSamples;
    using ghz::source::TeukolskyField;
    using ghz::source::TeukolskySourceCoefficients;
    using ghz::source::TeukolskySourceComponents;
    using ghz::source::TeukolskySourceOperators;
    using spectral::SpectralFieldVectorized;
    using teuk::Complex;
    using teuk::Real;

    struct TestFailure : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    bool approx_real(Real a, Real b, Real tol = Real(1e-11)) {
        return teuk::Fabs(a - b) <= tol;
    }

    bool approx_complex(const Complex& a, const Complex& b, Real tol = Real(1e-11)) {
        return approx_real(a.real(), b.real(), tol) &&
               approx_real(a.imag(), b.imag(), tol);
    }

    void expect_true(bool cond, const std::string& where) {
        if (!cond) {
            throw TestFailure(where);
        }
    }

    void expect_complex_eq(const Complex& got,
                           const Complex& expected,
                           const std::string& where,
                           Real tol = Real(1e-11)) {
        if (!approx_complex(got, expected, tol)) {
            std::ostringstream os;
            os << where << ": got " << got << ", expected " << expected;
            throw TestFailure(os.str());
        }
    }

    TeukolskyField make_field(size_t Nr,
                              size_t Nz,
                              int p,
                              int q,
                              Complex base) {
        SpectralFieldVectorized<ghp::GHPScalar<Complex>>::Modes modes{2, 0, 0};
        TeukolskyField field(Nr,
                             Nz,
                             modes,
                             ghp::GHPScalar<Complex>(teuk::zeroC, p, q),
                             p,
                             q);
        field.set_omega_mk(Real(3.25));

        for (size_t ir = 0; ir < Nr; ++ir) {
            for (size_t iz = 0; iz < Nz; ++iz) {
                const Complex offset(Real(ir + 1), Real(2 * iz + 1));
                field(ir, iz) = ghp::GHPScalar<Complex>(base + offset, p, q);
            }
        }

        return field;
    }

    TeukolskyField make_constant_coeff(size_t Nr,
                                       size_t Nz,
                                       int p,
                                       int q,
                                       Complex value) {
        SpectralFieldVectorized<ghp::GHPScalar<Complex>>::Modes modes{0, 0, 0};
        TeukolskyField field(Nr,
                             Nz,
                             modes,
                             ghp::GHPScalar<Complex>(teuk::zeroC, p, q),
                             p,
                             q);

        for (size_t ir = 0; ir < Nr; ++ir) {
            for (size_t iz = 0; iz < Nz; ++iz) {
                field(ir, iz) = ghp::GHPScalar<Complex>(value, p, q);
            }
        }

        return field;
    }

    TeukolskySourceOperators make_constant_operators(Complex thorn_factor,
                                                     Complex eth_factor) {
        TeukolskySourceOperators ops;
        ops.thorn = [thorn_factor](const TeukolskyField& in,
                                   TeukolskyField& out) {
            for (size_t ir = 0; ir < in.Nr(); ++ir) {
                for (size_t iz = 0; iz < in.Nz(); ++iz) {
                    out(ir, iz) = ghp::GHPScalar<Complex>(
                            thorn_factor * in(ir, iz).value(),
                            in.p() + 1,
                            in.q() + 1);
                }
            }
        };
        ops.eth = [eth_factor](const TeukolskyField& in,
                               TeukolskyField& out) {
            for (size_t ir = 0; ir < in.Nr(); ++ir) {
                for (size_t iz = 0; iz < in.Nz(); ++iz) {
                    out(ir, iz) = ghp::GHPScalar<Complex>(
                            eth_factor * in(ir, iz).value(),
                            in.p(),
                            in.q() - 2);
                }
            }
        };
        return ops;
    }

    void test_s0_teukolsky_formula() {
        const size_t Nr = 2;
        const size_t Nz = 3;

        const Complex thorn_factor(Real(2), Real(0.25));
        const Complex eth_factor(Real(-1), Real(0.5));

        const Complex rho(Real(0.2), Real(-0.1));
        const Complex rho_bar(Real(0.2), Real(0.1));
        const Complex tau(Real(-0.3), Real(0.4));
        const Complex taup_bar(Real(0.5), Real(-0.2));

        TeukolskySourceComponents components{
                make_field(Nr, Nz, 1, 3, Complex(Real(1), Real(0))),
                make_field(Nr, Nz, 0, 0, Complex(Real(2), Real(-1))),
                make_field(Nr, Nz, -1, -3, Complex(Real(-1), Real(2)))};

        TeukolskySourceCoefficients coeffs{
                make_constant_coeff(Nr, Nz, 1, 1, rho),
                make_constant_coeff(Nr, Nz, 1, 1, rho_bar),
                make_constant_coeff(Nr, Nz, 0, -2, tau),
                make_constant_coeff(Nr, Nz, 0, -2, taup_bar)};

        const TeukolskyField out = ghz::source::apply_s0_teukolsky_source(
                components,
                coeffs,
                make_constant_operators(thorn_factor, eth_factor));

        expect_true(out.p() == 1 && out.q() == -1,
                    "test_s0_teukolsky_formula: wrong output weight");
        expect_true(out.m() == 2 && out.has_omega_mk(),
                    "test_s0_teukolsky_formula: mode metadata not preserved");

        for (size_t ir = 0; ir < Nr; ++ir) {
            for (size_t iz = 0; iz < Nz; ++iz) {
                const Complex Tll = components.Tll(ir, iz).value();
                const Complex Tlm = components.Tlm_sym(ir, iz).value();
                const Complex Tmm = components.Tmm(ir, iz).value();

                const Complex A = (thorn_factor - Real(2) * rho_bar) * Tlm -
                                  (eth_factor - taup_bar) * Tll;
                const Complex B = (eth_factor - Real(2) * taup_bar) * Tlm -
                                  (thorn_factor - rho_bar) * Tmm;
                const Complex expected =
                        Real(0.5) *
                        ((eth_factor - taup_bar - Real(4) * tau) * A +
                         (thorn_factor - Real(4) * rho - rho_bar) * B);

                expect_complex_eq(out(ir, iz).value(),
                                  expected,
                                  "test_s0_teukolsky_formula: value mismatch");
            }
        }
    }

    void test_angular_projection() {
        SpectralFieldVectorized<ghp::GHPScalar<Complex>>::Modes modes{3, 0, 0};
        TeukolskyField source(2,
                              2,
                              modes,
                              ghp::GHPScalar<Complex>(teuk::zeroC, 1, -1),
                              1,
                              -1);

        source(0, 0) = ghp::GHPScalar<Complex>(Complex(Real(1), Real(2)), 1, -1);
        source(0, 1) = ghp::GHPScalar<Complex>(Complex(Real(3), Real(4)), 1, -1);
        source(1, 0) = ghp::GHPScalar<Complex>(Complex(Real(-1), Real(1)), 1, -1);
        source(1, 1) = ghp::GHPScalar<Complex>(Complex(Real(2), Real(-3)), 1, -1);

        AngularModeSamples harmonic;
        harmonic.spin = 2;
        harmonic.ell = 4;
        harmonic.m = 3;
        harmonic.values = {Complex(Real(1), Real(0)),
                           Complex(Real(0), Real(1))};

        const std::vector<Real> weights{Real(0.25), Real(0.75)};
        const std::vector<Complex> projected =
                ghz::source::project_angular_mode(source, harmonic, weights);

        const Complex row0 = teuk::twoPi *
                             (weights[0] * source(0, 0).value() *
                                      std::conj(harmonic.values[0]) +
                              weights[1] * source(0, 1).value() *
                                      std::conj(harmonic.values[1]));
        const Complex row1 = teuk::twoPi *
                             (weights[0] * source(1, 0).value() *
                                      std::conj(harmonic.values[0]) +
                              weights[1] * source(1, 1).value() *
                                      std::conj(harmonic.values[1]));

        expect_complex_eq(projected[0], row0, "test_angular_projection: row 0");
        expect_complex_eq(projected[1], row1, "test_angular_projection: row 1");
    }

    void test_bl_kinnersley_coefficients_feed_teukolsky_source() {
        KerrParams params(Real(1), Real(0.7));
        KerrMetric metric(params);
        CoordinateHelper coords(metric, CoordType::BoyerLindquist);
        KinnersleyTetrad<BLCoords> tetrad(metric, coords);

        const std::vector<Real> r_nodes{Real(4.0), Real(6.0)};
        const std::vector<Real> theta_nodes{
                Real(1.2),
                Real(2.0)};
        const std::vector<Real> z_nodes{
                teuk::Cos(theta_nodes[0]),
                teuk::Cos(theta_nodes[1])};

        const auto coeffs =
                ghz::source::make_bl_teukolsky_source_coefficients_from_theta(
                        tetrad,
                        r_nodes,
                        theta_nodes);

        const Real a = params.a;
        const Real sqrt2 = teuk::Sqrt(Real(2));

        for (size_t ir = 0; ir < r_nodes.size(); ++ir) {
            for (size_t iz = 0; iz < z_nodes.size(); ++iz) {
                const Real r = r_nodes[ir];
                const Real z = z_nodes[iz];
                const Real s1 = teuk::Sqrt(Real(1) - z * z);
                const Real sigma = r * r + a * a * z * z;

                const Complex rho = -Real(1) / (r - teuk::I * a * z);
                const Complex pi = teuk::I * a * s1 * rho * rho / sqrt2;
                const Complex tau = -teuk::I * a * s1 * rho * std::conj(rho) /
                                    sqrt2;
                const Complex taup_bar = std::conj(-pi);

                expect_complex_eq(coeffs.rho(ir, iz).value(),
                                  rho,
                                  "test_bl_kinnersley_coefficients: rho");
                expect_complex_eq(coeffs.rho_bar(ir, iz).value(),
                                  std::conj(rho),
                                  "test_bl_kinnersley_coefficients: rho_bar");
                expect_complex_eq(coeffs.tau(ir, iz).value(),
                                  tau,
                                  "test_bl_kinnersley_coefficients: tau");
                expect_complex_eq(coeffs.tau_prime_bar(ir, iz).value(),
                                  taup_bar,
                                  "test_bl_kinnersley_coefficients: tau_prime_bar");
            }
        }
    }

} // namespace

int main() {
    try {
        test_s0_teukolsky_formula();
        test_angular_projection();
        test_bl_kinnersley_coefficients_feed_teukolsky_source();
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }

    std::cout << "teukolsky_source_unit_tests: all tests passed\n";
    return 0;
}
