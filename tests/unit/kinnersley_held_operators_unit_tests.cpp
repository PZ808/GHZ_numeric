//
// Created by Peter Zimmerman on 07.03.26.
//

/**
 * @file test_kinnersley_held_operators.cpp
 *
 *  @brief Unit tests for Kinnersley/Held operators on spectral fields. \n
 *  @details
 *  This file is intended to test: \n
 *   1. basic action of edthH / edthBarH / thornPH on a single RSlice  \n
 *   2. consistency between dz_Dmatrix and dz_barycentric backends  \n
 *   3. spin-weight shifts n
 *   4. commutators of operators  \n
 *
 *
 */

#include "ghz/core/GhzTypes.hpp"
#include "ghz/geom/KerrParams.hpp"
#include "ghz/geom/KerrMetric.hpp"
#include "ghz/geom/Coords.hpp"
#include "ghz/geom/KerrCharts.hpp"
#include "ghz/spectral/SpectralCoordinateMaps.hpp"

#include "ghz/spectral/KinnersleySpectralHeldOperators.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/ghp/GHPScalars.hpp"

#include <cmath>
#include <complex>
#include <exception>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

    using teuk::Real;
    using teuk::Complex;
    using spectral::SpectralDiffer;
    using spectral::SpectralGHPVectorized;
    using ghp::HeldBackgroundFieldsVectorized;

    ghz::numeric::Domain domain{20, 100.0};

    template <typename T>
    using GHP = GHPScalar<T>;

    constexpr Real mass = 1.0;
    constexpr Real spin = 0.5;
    const int p_test = 2, q_test = 2;
    const int m_test = 2, kr_test = 0, kz_test = 0;
    const int slice_idx = 3;
    constexpr size_t Nz = 33;
    constexpr size_t Nr = 33;


    struct TestFailure : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    inline void require_true(bool cond, const std::string& msg) {
        if (!cond) throw TestFailure(msg);
    }

    inline Real abs_err(const Complex& a, const Complex& b) {
        return std::abs(a - b);
    }

    inline void require_near(const Complex& got, const Complex& expected,
                             Real tol, const std::string& msg) {
        const Real err = abs_err(got, expected);
        if (err > tol) {
            std::ostringstream oss;
            oss << msg << "\n"
                << "  got      = " << got << "\n"
                << "  expected = " << expected << "\n"
                << "  |err|    = " << err << "\n"
                << "  tol      = " << tol;
            throw TestFailure(oss.str());
        }
    }

    template <typename Field>
    Real max_abs_value_field(const Field& f) {
        Real emax = Real(0);
        for (size_t ir = 0; ir < f.Nr(); ++ir) {
            for (size_t iz = 0; iz < f.Nz(); ++iz) {
                emax = std::max(emax, std::abs(f(ir, iz).value()));
            }
        }
        return emax;
    }

    inline void require_near_real(Real got, Real expected,
                                  Real tol, const std::string& msg) {
        const Real err = std::abs(got - expected);
        if (err > tol) {
            std::ostringstream oss;
            oss << msg << "\n"
                << "  got      = " << got << "\n"
                << "  expected = " << expected << "\n"
                << "  |err|    = " << err << "\n"
                << "  tol      = " << tol;
            throw TestFailure(oss.str());
        }
    }

    template <typename SliceA, typename SliceB>
    Real max_abs_diff(const SliceA& a, const SliceB& b) {
        require_true(a.size() == b.size(), "max_abs_diff: slice size mismatch");
        Real emax = Real(0);
        for (size_t i = 0; i < a.size(); ++i) {
            emax = std::max(emax, std::abs(a[i].value() - b[i].value()));
        }
        return emax;
    }

    template <typename Slice>
    void require_uniform_pq(const Slice& s, int p, int q, const std::string& label) {
        for (size_t i = 0; i < s.size(); ++i) {
            if (s[i].p() != p || s[i].q() != q) {
                std::ostringstream oss;
                oss << label << ": wrong (p,q) at i=" << i
                    << ", got (" << s[i].p() << "," << s[i].q() << ")"
                    << ", expected (" << p << "," << q << ")";
                throw TestFailure(oss.str());
            }
        }
    }

    template <typename Slice>
    void zero_slice(Slice& s, int p, int q) {
        for (size_t i = 0; i < s.size(); ++i) {
            s[i] = GHP<Complex>(Complex(0.0, 0.0), p, q);
        }
    }

    template <typename Slice>
    void copy_slice(const Slice& in, Slice& out) {
        require_true(in.size() == out.size(), "copy_slice: size mismatch");
        for (size_t i = 0; i < in.size(); ++i) out[i] = in[i];
    }

    //inline KinnersleyTetrad<OutgoingCoords>
    //build_test_tetrad(KerrMetric kerr, CoordinateHelper coords) {
     //   return KinnersleyTetrad<OutgoingCoords>(kerr, coords);
   // }

    inline OutgoingCoords build_test_Xout() {
        return OutgoingCoords(10.0, Real(12345), Real(.23456), M_PI_2);
    }

    inline ghp::GHPBackgroundFieldsVectorized
    build_test_ghp_background(const SpectralDiffer& diff,
                              KinnersleyTetrad<OutgoingCoords>& tetrad,
                              OutgoingCoords& Xout, const ghz::numeric::AffineMap1D& r_map) {

        auto physical_nodes = r_map.toPhysicalNodes(diff.cl_nodes());
        return ghp::build_ghp_fields_vectorized( tetrad, physical_nodes, diff.lgl_nodes(), Xout);
    }

    inline ghp::HeldBackgroundFieldsVectorized<OutgoingCoords>
    build_test_held_background(const SpectralDiffer& diff,
                               KinnersleyTetrad<OutgoingCoords>& tetrad,
                               OutgoingCoords& Xout) {
        return ghp::build_held_fields_vectorized( tetrad, diff.lgl_nodes(), Xout).fields;
    }

    struct TestObjects {
        KerrParams params{mass, spin};
        KerrMetric kerr{params};
        CoordinateHelper coords{kerr};
        SpectralDiffer diff;
        KinnersleyTetrad<OutgoingCoords> tetrad;
        OutgoingCoords Xout;
        ghz::numeric::AffineMap1D r_map;
        ghp::GHPBackgroundFieldsVectorized bkg_ghp_fields;
        ghp::HeldBackgroundFieldsVectorized<OutgoingCoords> bkg_held_fields;
        KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer> ops;

        TestObjects(size_t nz, size_t nr)
                : diff(nz, nr),
                  tetrad(kerr, coords),
                  Xout(build_test_Xout()),
                  r_map(domain),
                  bkg_ghp_fields(build_test_ghp_background(diff, tetrad, Xout, r_map)),
                  bkg_held_fields(build_test_held_background(diff, tetrad, Xout)),
                  ops(diff, bkg_held_fields, tetrad, r_map)
        {}
    };

// Build a 2D field and return it.
// Here initialize:
//   - dimensions Nz, Nr
//   - mode metadata
//   - omega if needed
//   - values f(r_i, z_j)
    SpectralGHPVectorized
    build_test_field(const SpectralDiffer& diff, const ghz::numeric::AffineMap1D& r_map,
                     int p, int q,
                     int m, int kr, int kz,
                     Real omega_mk)
    {
        using Scalar = GHPScalar<Complex>;

        ghp::GHPFieldVectorized f_test(
                diff.Nr(), diff.Nz(),
                Scalar(Complex(0.0, 0.0), p, q),
                p, q
        );

        auto r_nodes = r_map.toPhysicalNodes(diff.cl_nodes());
        // choose a function with some r and z dependence and a nontrivial complex phase to test the operators on
        auto fillerInsane = [&](Real r, Real z) -> Complex {
            return  std::pow(1.0 - z*z, 3/2)*(Real(1.0) + Real(0.2)*r + Real(0.05)*r*r+Real(2.3)*r*r*r*r)
                   * (Real(1.0) + Real(0.3)*z + Real(0.2)*z*z+Real(10)*z*z*z)
                   * std::exp(Complex(0.0, Real(0.4)*r*r - Real(0.3)*z));
        };

        auto fillerDumb = [&](Real r, Real z) -> Complex {
            return (Real(1.0) + Real(0.01)*r + Real(0.0002)*r*r)
                   * (Real(1.0) + Real(0.3)*z + Real(0.2)*z*z)
                   * std::exp(Complex(0.0, Real(0.01)*r - Real(0.1)*z));
        };

        auto filler = [&](Real r, Real z) -> Complex {
            return (Real(1.0) + Real(0.01)*r + Real(0.0002)*r*r)
                   * (Real(1.0) + Real(0.3)*z + Real(0.2)*z*z);
        };
        auto fillerConst = [&](Real r, Real z) -> Complex { return Complex(1,0); };

        f_test.fill_nodes(fillerConst, r_nodes, diff.lgl_nodes());

        SpectralGHPVectorized::Modes modes{m, kr, kz};

        SpectralGHPVectorized out(diff.Nr(), diff.Nz(), modes, f_test);
        out.set_omega_mk(omega_mk);

        return out;
    }

    void debug_edthH_on_constant()
    {
        TestObjects T(Nz, Nr);

        auto field = build_test_field(
                T.diff, T.r_map,
                /*p=*/0, /*q=*/0,
                /*m=*/0, /*kr=*/0, /*kz=*/0,
                /*omega=*/Real(0));

        // overwrite with literal constant 1
        for (size_t ir = 0; ir < field.Nr(); ++ir) {
            for (size_t iz = 1; iz < field.Nz()-1; ++iz) {
                field(ir, iz) = GHP<Complex>(Complex(1.0, 0.0), 0, 0);
            }
        }

        auto eth_f = SpectralGHPVectorized(
                field.Nr(), field.Nz(), field.modes(),
                GHP<Complex>(Complex(0.0, 0.0), 0, -2), 0, -2);

        auto ethbar_f = SpectralGHPVectorized(
                field.Nr(), field.Nz(), field.modes(),
                GHP<Complex>(Complex(0.0, 0.0), -2, 0), -2, 0);

        auto ethbar_eth_f = SpectralGHPVectorized(
                field.Nr(), field.Nz(), field.modes(),
                GHP<Complex>(Complex(0.0, 0.0), -2, -2), -2, -2);

        auto eth_ethbar_f = SpectralGHPVectorized(
                field.Nr(), field.Nz(), field.modes(),
                GHP<Complex>(Complex(0.0, 0.0), -2, -2), -2, -2);

        if (field.has_omega_mk()) {
            eth_f.set_omega_mk(field.omega_mk());
            ethbar_f.set_omega_mk(field.omega_mk());
            ethbar_eth_f.set_omega_mk(field.omega_mk());
            eth_ethbar_f.set_omega_mk(field.omega_mk());
        }

        T.ops.edthH_inplace(field, eth_f);
        T.ops.edthBarH_inplace(field, ethbar_f);

        T.ops.edthBarH_inplace(eth_f, ethbar_eth_f);
        T.ops.edthH_inplace(ethbar_f, eth_ethbar_f);

        std::cout << "\n[DEBUG constant mode]\n";
        std::cout << "max |edthH(f)|              = " << max_abs_value_field(eth_f) << "\n";
        std::cout << "max |edthBarH(f)|           = " << max_abs_value_field(ethbar_f) << "\n";
        std::cout << "max |edthBarH(edthH(f))|    = " << max_abs_value_field(ethbar_eth_f) << "\n";
        std::cout << "max |edthH(edthBarH(f))|    = " << max_abs_value_field(eth_ethbar_f) << "\n";

        // Also print a few sample points
        for (size_t ir = 0; ir < std::min<size_t>(field.Nr(), 3); ++ir) {
            for (size_t iz = 0; iz < std::min<size_t>(field.Nz(), 3); ++iz) {
                std::cout << "ir=" << ir << ", iz=" << iz
                          << "  edth=" << eth_f(ir, iz).value()
                          << "  edthBar=" << ethbar_f(ir, iz).value()
                          << "  ethbar_eth=" << ethbar_eth_f(ir, iz).value()
                          << "  eth_ethbar=" << eth_ethbar_f(ir, iz).value()
                          << "\n";
            }
        }
    }

// -----------------------------------------------------------------------------
// Helpers to apply operators to a single RSlice
// -----------------------------------------------------------------------------
    SpectralGHPVectorized::RSlice make_output_rslice_like(
            const SpectralGHPVectorized::RSlice& in,
            std::vector<GHP<Complex>>& backing,
            int p, int q)
    {
        backing.assign(in.size(), GHP<Complex>(Complex(0.0, 0.0), p, q));

        return SpectralGHPVectorized::RSlice(
                backing.data(),
                in.size(),
                in.modes_,
                in.r_index(),
                in.has_omega_mk() ? in.omega_mk() : Real(0),
                in.has_omega_mk()
        );
    }

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

    void test_edth_spin_shift() {

        TestObjects T(Nz, Nr);

        auto field = build_test_field(T.diff, T.r_map, p_test, q_test,
                                      m_test, kr_test, kz_test, Real(0.5));
        std::vector<GHP<Complex>> out_buf(field.slice_R(slice_idx).size());
        auto in = field.slice_R(slice_idx);
        auto out = make_output_rslice_like(in, out_buf, 0, 0);

        T.ops.edthH_inplace_RSliceV(in, out);

        // Eth: (p,q) -> (p, q-2)
        require_uniform_pq(out, p_test, q_test-2, "edthH spin shift");
    }

    void test_edthbar_spin_shift() {

        TestObjects T(Nz, Nr);

        auto field = build_test_field(T.diff, T.r_map, p_test, q_test,
                                      m_test, kr_test, kz_test, Real(0.5));
        auto in = field.slice_R(slice_idx);
        std::vector<GHP<Complex>> out_buf(in.size());
        auto out = make_output_rslice_like(in, out_buf, 0, 0);

        T.ops.edthBarH_inplace_RSliceV(in, out);

        // EthBar: (p,q) -> (p-2, q)
        require_uniform_pq(out, p_test-2, q_test, "edthBarH spin shift");
    }

    void test_thornPH_spin_shift() {

        TestObjects T(Nz, Nr);

        auto field = build_test_field(T.diff, T.r_map, p_test, q_test,
                                      m_test, kr_test, kz_test, Real(0.5));
        auto in = field.slice_R(slice_idx);
        std::vector<GHP<Complex>> out_buf(in.size());
        auto out = make_output_rslice_like(in, out_buf, 0, 0);

        T.ops.thornPH_inplace_RSliceV(in, out);

        // thornPH: (p,q) -> (p-1, q-1)
        require_uniform_pq(out, p_test-1, q_test-1, "thornPH spin shift");
    }

    void test_edth_dmatrix_vs_bary() {

        TestObjects T(Nz, Nr);

        auto field = build_test_field(T.diff, T.r_map, p_test, q_test,
                                      m_test, kr_test, kz_test, Real(0.5));

        auto in = field.slice_R(slice_idx);

        std::vector<GHP<Complex>> out1_buf(in.size());
        std::vector<GHP<Complex>> out2_buf(in.size());

        auto out1 = make_output_rslice_like(in, out1_buf, 0, 0);
        auto out2 = make_output_rslice_like(in, out2_buf, 0, 0);

        T.ops.edthH_dmat_inplace_RSliceV(in, out1);
        T.ops.edthH_bary_inplace_RSliceV(in, out2);

        const Real err = max_abs_diff(out1, out2);
        require_near_real(err, Real(0), Real(1e-8), "edth D-matrix vs barycentric");
    }

    void test_thornPH_against_expected_frequency_factor() {

        const Real omega = Real(1.234);
        TestObjects T(Nz, Nr);
        auto field = build_test_field(T.diff, T.r_map, p_test, q_test,
                                      m_test, kr_test, kz_test,omega);
        auto in = field.slice_R(slice_idx);

        std::vector<GHP<Complex>> out_buf(in.size());
        auto out = make_output_rslice_like(in, out_buf, 0, 0);

        T.ops.thornPH_inplace_RSliceV(in, out);

        const Complex minus_i_omega = -Complex(0.0, 1.0) * omega;

        for (size_t i = 0; i < in.size(); ++i) {
            require_near(out[i].value(),
                         minus_i_omega * in[i].value(),
                         Real(1e-13),
                         "thornPH frequency factor at i=" + std::to_string(i));
        }

        require_uniform_pq(out, p_test - 1, q_test - 1, "thornPH expected spin shift");
    }

// -----------------------------------------------------------------------------
// Commutator scaffolding
// -----------------------------------------------------------------------------

    template <typename ApplyA, typename ApplyB>
    Real commutator_norm(const SpectralGHPVectorized::RSlice& in,
                         ApplyA A,
                         ApplyB B)
    {
        const size_t N = in.size();

        std::vector<GHP<Complex>> a_buf(N);
        std::vector<GHP<Complex>> b_buf(N);
        std::vector<GHP<Complex>> ab_buf(N);
        std::vector<GHP<Complex>> ba_buf(N);

        // We do not know output (p,q) in advance for arbitrary A,B,
        // so initialize with input metadata and let operators overwrite.
        SpectralGHPVectorized::RSlice a(
                a_buf.data(), N, in.modes_, in.r_index(),
                in.has_omega_mk() ? in.omega_mk() : Real(0),
                in.has_omega_mk()
        );

        SpectralGHPVectorized::RSlice b(
                b_buf.data(), N, in.modes_, in.r_index(),
                in.has_omega_mk() ? in.omega_mk() : Real(0),
                in.has_omega_mk()
        );

        SpectralGHPVectorized::RSlice ab(
                ab_buf.data(), N, in.modes_, in.r_index(),
                in.has_omega_mk() ? in.omega_mk() : Real(0),
                in.has_omega_mk()
        );

        SpectralGHPVectorized::RSlice ba(
                ba_buf.data(), N, in.modes_, in.r_index(),
                in.has_omega_mk() ? in.omega_mk() : Real(0),
                in.has_omega_mk()
        );

        for (size_t i = 0; i < N; ++i) {
            a[i]  = GHP<Complex>(Complex(0.0, 0.0), in[i].p(), in[i].q());
            b[i]  = GHP<Complex>(Complex(0.0, 0.0), in[i].p(), in[i].q());
            ab[i] = GHP<Complex>(Complex(0.0, 0.0), in[i].p(), in[i].q());
            ba[i] = GHP<Complex>(Complex(0.0, 0.0), in[i].p(), in[i].q());
        }

        A(in, a);
        B(a, ab);

        B(in, b);
        A(b, ba);

        return max_abs_diff(ab, ba);
    }
    template <typename ApplyA, typename ApplyB, typename ApplyRHS>
    Real commutator_rhs_norm_field(const SpectralGHPVectorized& in,
                                   ApplyA A,
                                   ApplyB B,
                                   ApplyRHS RHS,
                                   int p_mid_A, int q_mid_A,
                                   int p_mid_B, int q_mid_B,
                                   int p_out, int q_out)
    {
        auto make_like = [&](int p, int q) {
            SpectralGHPVectorized out(
                    in.Nr(), in.Nz(),
                    in.modes(),
                    GHP<Complex>(Complex(0.0, 0.0), p, q),
                    p, q
            );
            if (in.has_omega_mk()) out.set_omega_mk(in.omega_mk());
            return out;
        };

        auto a   = make_like(p_mid_A, q_mid_A);
        auto b   = make_like(p_mid_B, q_mid_B);
        auto ab  = make_like(p_out, q_out);
        auto ba  = make_like(p_out, q_out);
        auto rhs = make_like(p_out, q_out);

        A(in, a);
        B(a, ab);

        B(in, b);
        A(b, ba);

        RHS(in, rhs);

        Real emax = Real(0);
        for (size_t ir = 0; ir < in.Nr(); ++ir) {
            for (size_t iz = 0; iz < in.Nz(); ++iz) {
                const Complex diff = ab(ir, iz).value() - ba(ir, iz).value() - rhs(ir, iz).value();
                emax = std::max(emax, std::abs(diff));
            }
        }
        return emax;
    }




    Real commutator_edthH_edthbarH_rhs_norm(const TestObjects& T,
                                          const SpectralGHPVectorized& f)
    {
        const int p = f.p();
        const int q = f.q();

        // [edth, edthbar] lands at (p-2, q-2)
        auto eth_ethbar = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-2, q-2), p-2, q-2
        );
        auto ethbar_eth = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-2, q-2), p-2, q-2
        );

        auto tmp_eth = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p, q-2), p, q-2
        );
        auto tmp_ethbar = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-2, q), p-2, q
        );

        // thorn(f) -> (p+1,q+1)
        auto thorn_f = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p+1, q+1), p+1, q+1
        );

        // thornPH(f) -> (p-1,q-1)
        auto thornP_f = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-1, q-1), p-1, q-1
        );

        auto rhs = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-2, q-2), p-2, q-2
        );

        if (f.has_omega_mk()) {
            tmp_eth.set_omega_mk(f.omega_mk());
            tmp_ethbar.set_omega_mk(f.omega_mk());
            thorn_f.set_omega_mk(f.omega_mk());
            thornP_f.set_omega_mk(f.omega_mk());
            eth_ethbar.set_omega_mk(f.omega_mk());
            ethbar_eth.set_omega_mk(f.omega_mk());
            rhs.set_omega_mk(f.omega_mk());
        }

        T.ops.edthH_inplace(f, tmp_eth);
        T.ops.edthBarH_inplace(tmp_eth, ethbar_eth);

        T.ops.edthBarH_inplace(f, tmp_ethbar);
        T.ops.edthH_inplace(tmp_ethbar, eth_ethbar);

        T.ops.thorn_inplace(f, thorn_f);
        T.ops.thornPHr_inplace(T.tetrad, f, thornP_f);

        const auto& held = T.bkg_held_fields;
        const auto& ghp = T.bkg_ghp_fields;

        Real emax = Real(0);
        for (size_t iz = 1; iz < f.Nz()-1; ++iz) {
            // Replace these with your exact stored Held fields / conventions
            const Complex Om0= held.OmH(iz).value();
            const Complex tau0 = held.tauH(iz).value();
            const Complex taub0 = held.tauH_bar(iz).value();
            const Complex rhop0 = held.rhopH(iz).value();
            const Complex rhopb0 = held.rhopH_bar(iz).value();
            const Complex psi0 = held.PsiH(iz).value();
            const Complex psib0 = held.PsiH_bar(iz).value();


            for (size_t ir = 0; ir < f.Nr(); ++ir) {
                const Complex rho = ghp.rho(ir, iz).value();
                const Complex rhob = std::conj(rho);
                const Complex rhop = ghp.rhop(ir, iz).value();
                const Complex rhopb = std::conj(ghp.rhop(ir, iz).value());
                //const Complex EdthH_taubar0 = -half*Om0*(rhop0+rhopb0) - half*(psi0-psib0);
                const Complex EdthH_taubar0 =
                        -Om0 * rhop0 + Real(0.5) * (psib0 - psi0);
                const Complex EdthH_rho = sqr(rho)*tau0;
                const Complex psi2 = cube(rho)*psi0;
                const Complex psi2b = std::conj(cube(rho)*psi0);
                // Held Eq. (4.8)
                const Complex Sig0 = rhop/rhob + psi2 / (Real(2)*rho) * (Real(1)/rho + Real(1)/rhob) + EdthH_taubar0 * rho + taub0 * EdthH_rho;
                const Complex Sig0b = std::conj(Sig0);

                const Complex thornCoeff = (rhopb-rhop)/(rho*rhob);
                const Complex thornpHCoeff = Om0;

                const Complex Pterm = rhop/rhob+psi2/(Real(2)*rho)*(Real(1)/rho+Real(1)/rhob);
                const Complex Qterm = rhopb/rho+psi2b/(Real(2)*rhob)*(Real(1)/rho+Real(1)/rhob);
                // Held Eq. (4.7)
                const Complex rhs_val =
                        thornCoeff * thorn_f(ir, iz).value()
                        + thornpHCoeff * thornP_f(ir, iz).value()
                        + (Real(p)*Pterm + Real(q)*Qterm) * f(ir, iz).value();

                rhs(ir, iz) = GHP<Complex>(rhs_val, p-2, q-2);

                const Complex lhs_val =
                        eth_ethbar(ir, iz).value() - ethbar_eth(ir, iz).value();

                emax = std::max(emax, std::abs(lhs_val - rhs_val));

            }
        }

        Real emax2 = Real(0);
        size_t worst_ir = 0, worst_iz = 0;
        Complex worst_lhs = 0.0;
        Complex worst_rhs = 0.0;
        Complex worst_res = 0.0;

        for (size_t iz = 2; iz < f.Nz()-4; ++iz) {
            const Complex Om0 = held.OmH(iz).value();
            const Complex tau0 = held.tauH(iz).value();
            const Complex taub0 = held.tauH_bar(iz).value();
            const Complex rhop0 = held.rhopH(iz).value();
            const Complex rhopb0 = held.rhopH_bar(iz).value();
            const Complex psi0 = held.PsiH(iz).value();
            const Complex psib0 = held.PsiH_bar(iz).value();

            for (size_t ir = 0; ir < f.Nr(); ++ir) {
                const Complex rho = ghp.rho(ir, iz).value();
                const Complex rhob = std::conj(rho);
                const Complex rhop = ghp.rhop(ir, iz).value();
                const Complex rhopb = std::conj(ghp.rhop(ir, iz).value());

                const Complex EdthH_taubar0 =
                        -half*Om0*(rhop0+rhopb0) - half*(psi0-psib0);
                const Complex EdthH_rho = sqr(rho)*tau0;
                const Complex psi2 = cube(rho)*psi0;
                const Complex Sig0 = rhop/rhob
                                     + psi2 / (Real(2)*rho) * (Real(1)/rho + Real(1)/rhob)
                                     + EdthH_taubar0*rho + taub0*EdthH_rho;

                const Complex thornCoeff = (rhopb-rhop)/(rho*rhob);
                const Complex thornpHCoeff = Om0;

                const Complex rhs_val =
                        thornCoeff * thorn_f(ir, iz).value()
                        + thornpHCoeff * thornP_f(ir, iz).value()
                        + (Real(p)*Sig0 - Real(q)*std::conj(Sig0)) * f(ir, iz).value();

                const Complex lhs_val =
                        eth_ethbar(ir, iz).value() - ethbar_eth(ir, iz).value();

                const Complex res = lhs_val - rhs_val;
                const Real err = std::abs(res);

                if (err > emax2) {
                    emax2 = err;
                    worst_ir = ir;
                    worst_iz = iz;
                    worst_lhs = lhs_val;
                    worst_rhs = rhs_val;
                    worst_res = res;
                }
            }
        }

        auto r_nodes = T.r_map.toPhysicalNodes(T.diff.cl_nodes());

        std::cout << "\n[DEBUG commutator worst point]\n";
        std::cout << "worst ir = " << worst_ir
                  << ", worst iz = " << worst_iz << "\n";
        std::cout << "r = " << r_nodes[worst_ir]
                  << ", z = " << T.diff.lgl_nodes()[worst_iz] << "\n";
        std::cout << "lhs = " << worst_lhs << "\n";
        std::cout << "rhs = " << worst_rhs << "\n";
        std::cout << "res = " << worst_res << "\n";



        return emax2;
    }

    Real commutator_thorn_thornPH_rhs_norm(const TestObjects& T,
                                           const SpectralGHPVectorized& f)
    {
        const int p = f.p();
        const int q = f.q();

        // [thorn, thornPH] lands back at (p, q)
        auto thorn_thornP = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p, q), p, q
        );
        auto thornP_thorn = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p, q), p, q
        );

        // intermediates
        auto tmp_thorn = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p+1, q+1), p+1, q+1
        );
        auto tmp_thornP = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-1, q-1), p-1, q-1
        );

        // thorn(f) for RHS coefficient * thorn(f)
        auto thorn_f = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p+1, q+1), p+1, q+1
        );

        auto rhs = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p, q), p, q
        );

        if (f.has_omega_mk()) {
            tmp_thorn.set_omega_mk(f.omega_mk());
            tmp_thornP.set_omega_mk(f.omega_mk());
            thorn_f.set_omega_mk(f.omega_mk());
            thorn_thornP.set_omega_mk(f.omega_mk());
            thornP_thorn.set_omega_mk(f.omega_mk());
            rhs.set_omega_mk(f.omega_mk());
        }

        // thorn(thornPH(f))
        T.ops.thornPHr_inplace(T.tetrad, f, tmp_thornP);
        T.ops.thorn_inplace(tmp_thornP, thorn_thornP);

        // thornPH(thorn(f))
        T.ops.thorn_inplace(f, tmp_thorn);
        T.ops.thornPHr_inplace(T.tetrad, tmp_thorn, thornP_thorn);

        // thorn(f) for RHS
        T.ops.thorn_inplace(f, thorn_f);

        const auto& held = T.bkg_held_fields;
        const auto& ghp = T.bkg_ghp_fields;

        Real emax = Real(0);
        for (size_t iz = 0; iz < f.Nz(); ++iz) {
            const Complex tau0 = held.tauH(iz).value();
            const Complex taub0 = held.tauH_bar(iz).value();
            const Complex psi0 = held.PsiH(iz).value();
            const Complex psib0 = held.PsiH_bar(iz).value();

            for (size_t ir = 0; ir < f.Nr(); ++ir) {
                const Complex rho = ghp.rho(ir, iz).value();
                const Complex rhob = std::conj(rho);

                // full tau
                const Complex tau = tau0*rho*rhob;
                const Complex taub = std::conj(tau);

                // Type D relation: tau' = -rho^2 * taub0
                const Complex taup = -sqr(rho)*taub0;
                const Complex taupb = std::conj(taup);

                // Psi_2 = rho^3 Psi^o
                const Complex psi2 = cube(rho)*psi0;
                const Complex psi2b = std::conj(psi2);

                const Complex is_Zero = (taub*taupb*rho + tau*taup*rhob)/(rhob*rho);
                // std::cout << "ir=" << ir << ", iz=" << iz << ", is_Zero=" << is_Zero << "\n"; // Held Eq. (4.7a)

                const Complex coeff = is_Zero
                        - half*psi2/rho
                        - half*psi2b/rhob;
                const Complex rhs_val = coeff*thorn_f(ir, iz).value();

                rhs(ir, iz) = GHP<Complex>(rhs_val, p, q);

                const Complex lhs_val =
                        thorn_thornP(ir, iz).value() - thornP_thorn(ir, iz).value();

                emax = std::max(emax, std::abs(lhs_val-rhs_val));
            }
        }

        return emax;
    }

    Real commutator_thornPH_edthbarH_rhs_norm(const TestObjects& T,
                                           const SpectralGHPVectorized& f)
    {
        const int p = f.p();
        const int q = f.q();

        // [thorn, thornPH] lands back at (p, q)
        auto thornPH_edthbarH = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-3, q-1), -3, -1
        );
        auto edthbarH_thornPH = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-3, q-1), -3, -1
        );

        // intermediates
        auto tmp_thornPH = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-1, q-1), p-1, q-1
        );
        auto tmp_edthbarH = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-2, q), p-2, q
        );


        auto rhs = SpectralGHPVectorized(
                f.Nr(), f.Nz(), f.modes(),
                GHP<Complex>(Complex(0.0, 0.0), p-3, q-1), p-3, q-1
        );

        if (f.has_omega_mk()) {
            tmp_edthbarH.set_omega_mk(f.omega_mk());
            tmp_thornPH.set_omega_mk(f.omega_mk());
            thornPH_edthbarH.set_omega_mk(f.omega_mk());
            edthbarH_thornPH.set_omega_mk(f.omega_mk());
            rhs.set_omega_mk(f.omega_mk());
        }

        // thorn(thornPH(f))
        T.ops.thornPHr_inplace(T.tetrad, f, tmp_thornPH);
        T.ops.edthBarH_inplace(tmp_thornPH, edthbarH_thornPH);

        // thornPH(thorn(f))
        T.ops.edthBarH_inplace(f, tmp_edthbarH);
        T.ops.thornPHr_inplace(T.tetrad, tmp_edthbarH, thornPH_edthbarH);

        const auto& held = T.bkg_held_fields;
        const auto& ghp = T.bkg_ghp_fields;

        Real emax = Real(0);
        for (size_t iz = 0; iz < f.Nz(); ++iz) {
            const Complex tau0 = held.tauH(iz).value();
            const Complex taub0 = held.tauH_bar(iz).value();
            const Complex psi0 = held.PsiH(iz).value();
            const Complex psib0 = held.PsiH_bar(iz).value();

            for (size_t ir = 0; ir < f.Nr(); ++ir) {
                const Complex rho = ghp.rho(ir, iz).value();
                const Complex rhob = std::conj(rho);

                // full tau
                const Complex tau = tau0*rho*rhob;
                const Complex taub = std::conj(tau);

                // Type D relation: tau' = -rho^2 * taub0
                const Complex taup = -sqr(rho)*taub0;
                const Complex taupb = std::conj(taup);

                // Psi_2 = rho^3 Psi^o
                const Complex psi2 = cube(rho)*psi0;
                const Complex psi2b = std::conj(psi2);

                const Complex rhs_val = Complex(0,0);

                rhs(ir, iz) = GHP<Complex>(rhs_val, p, q);

                const Complex lhs_val =
                        thornPH_edthbarH(ir, iz).value() - edthbarH_thornPH(ir, iz).value();

                emax = std::max(emax, std::abs(lhs_val-rhs_val));
            }
        }

        return emax;
    }


    void test_commutator_edthH_edthbarH_rhs()
    {
        TestObjects T(Nz, Nr);

        auto field = build_test_field(T.diff, T.r_map,
                                      p_test, q_test,
                                      m_test, kr_test, kz_test,
                                      Real(0.8));



        const Real err = commutator_edthH_edthbarH_rhs_norm(T, field);
        std::cout << "[INFO] [edthH, edthBarH] residual = " << err << "\n";

        require_near_real(
                err,
                Real(0),
                Real(1e-8),
                "[edthH,edthBarH]f - (((rhopb-rhop)/(rho*rhob))*thorn(f) + Om0*thornPH(f) + (p*Sig0 - q*conj(Sig0))*f)"
        );
    }

    void test_commutator_thorn_thornPH_rhs()
    {
        TestObjects T(Nz, Nr);

        auto field = build_test_field(T.diff, T.r_map, p_test, q_test,
                                      m_test, kr_test, kz_test,
                                      Real(0.8));

        const Real err = commutator_thorn_thornPH_rhs_norm(T, field);

        std::cout << "[INFO] [thorn, thornPH] residual = " << err << "\n";
    }


    void test_commutator_thornPH_edthbarH_rhs()
    {
        TestObjects T(Nz, Nr);

        auto field = build_test_field(T.diff, T.r_map, p_test, q_test,
                                      m_test, kr_test, kz_test,
                                      Real(0.8));

        const Real err = commutator_thornPH_edthbarH_rhs_norm(T, field);

        std::cout << "[INFO] [thornPH, edthbarH] residual = " << err << "\n";
    }

    void test_commutators_reduced_rslices() {

        TestObjects T(Nz, Nr);

        auto field = build_test_field(
                T.diff, T.r_map, /*p=*/p_test, /*q=*/q_test,
                /*m=*/m_test, /*kr=*/kr_test, /*kz=*/kz_test,
                Real(0.8));
        auto in = field.slice_R(slice_idx);

        // -----------------------------------------------------------------
        // Commutators to test will depend on the identities
        // Example:
        //   [thornPH, edthH]
        // -----------------------------------------------------------------

        auto ThornPHeld = [&](const auto& src, auto& dst) {
            T.ops.thornPH_inplace_RSliceV(src, dst);
        };

        auto EdthHeld = [&](const auto& src, auto& dst) {
            T.ops.edthH_inplace_RSliceV(src, dst);
        };

        auto EdthBarHeld = [&](const auto& src, auto& dst) {
            T.ops.edthBarH_inplace_RSliceV(src, dst);
        };

        const Real comm_thornp_eth    = commutator_norm(in, ThornPHeld, EdthHeld);
        const Real comm_thornp_ethbar = commutator_norm(in, ThornPHeld, EdthBarHeld);

        // Replace tolerances/expected values once you decide the exact
        // commutator identities you want to enforce.
        //
        // If a commutator should vanish:
        // require_near_real(comm_thorn_eth, 0.0, 1e-10, "[thornPH,edthH]");
        //
        // If it should equal another operator, you’ll want a slightly different
        // helper comparing [A,B]f against Cf.

        std::cout << "[INFO] [thornPH, edthH] max-norm    = " << comm_thornp_eth << "\n";
        std::cout << "[INFO] [thornPH, edthBarH] max-norm = " << comm_thornp_ethbar << "\n";


        require_near_real(comm_thornp_eth,
                          Real(0),
                          Real(1e-7),
                          "[thornPH, edthH]");

        require_near_real(comm_thornp_ethbar,
                          Real(0),
                          Real(1e-7),
                          "[thornPH, edthBarH]");
    }

    using TestFn = void(*)();

    struct TestCase {
        const char* name;
        TestFn fn;
    };

} // namespace


int main() {
    const std::vector<TestCase> tests = {
            {"edth_spin_shift", test_edth_spin_shift},
            {"edthbar_spin_shift", test_edthbar_spin_shift},
            {"thornPH_spin_shift", test_thornPH_spin_shift},
            {"edth_dmatrix_vs_bary", test_edth_dmatrix_vs_bary},
            {"thornPH_expected_frequency_factor", test_thornPH_against_expected_frequency_factor},
            {"commutators_rslices", test_commutators_reduced_rslices},
            {"commutator_thorn_thornPH_rhs", test_commutator_thorn_thornPH_rhs},
            {"commutator_edth_edthbar_rhs", test_commutator_edthH_edthbarH_rhs},
            {"commutator_thornPH_edthbar_rhs", test_commutator_thornPH_edthbarH_rhs}
    };

    int passed = 0;
    for (const auto& t : tests) {
        try {
            t.fn();
            std::cout << "[PASS] " << t.name << "\n";
            ++passed;
        } catch (const std::exception& e) {
            std::cerr << "[FAIL] " << t.name << "\n" << e.what() << "\n";
            return EXIT_FAILURE;
        } catch (...) {
            std::cerr << "[FAIL] " << t.name << "\nUnknown exception\n";
            return EXIT_FAILURE;
        }
    }

    std::cout << "\nAll " << passed << " Kinnersley/Held tests passed.\n";
    return EXIT_SUCCESS;
}
