//
// effective_source_unit_test.cpp
//
// Plain unit tests for ghz::source::EffectiveSourceArchive / EffectiveSourceSampler
//

#include "ghz/source/EffectiveSource.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

    using ghz::source::Axis;
    using ghz::source::Component;
    using ghz::source::EffectiveSourceArchive;
    using ghz::source::EffectiveSourceSampler;
    using ghz::source::ModeData;
    using ghz::source::PatchSide;
    using ghz::source::SourceMetadata;
    using ghz::source::YStorage;

    using teuk::Real;
    using teuk::Complex;

    struct TestFailure : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    template <typename T>
    void expect_true(bool cond, const T& msg) {
        if (!cond) throw TestFailure(msg);
    }

    bool approx_real(Real a, Real b, Real tol = Real(1e-12)) {
        return teuk::Fabs(a - b) <= tol;
    }

    bool approx_complex(const Complex& a, const Complex& b, Real tol = Real(1e-12)) {
        return approx_real(a.real(), b.real(), tol) &&
               approx_real(a.imag(), b.imag(), tol);
    }

    void expect_real_eq(Real got, Real expected, const std::string& where,
                        Real tol = Real(1e-12)) {
        if (!approx_real(got, expected, tol)) {
            std::ostringstream os;
            os << where << ": got " << got << ", expected " << expected;
            throw TestFailure(os.str());
        }
    }

    void expect_complex_eq(const Complex& got, const Complex& expected,
                           const std::string& where,
                           Real tol = Real(1e-12)) {
        if (!approx_complex(got, expected, tol)) {
            std::ostringstream os;
            os << where << ": got " << got << ", expected " << expected;
            throw TestFailure(os.str());
        }
    }

    // -------------------------------------------------------------------------
    // Fake archive
    //
    // Raw storage:
    //   Left patch  : X = [-2, -1]
    //   Right patch : X = [ 1,  2]
    //   Y = [0, 1, 2]   (half-grid, includes Y=0)
    //
    // For component ll only:
    //   Re f_m(X,Y) = 100 m + 10 X + Y
    //   Im f_m(X,Y) = 1000 m + 100 X + 10 Y
    //
    // All other components are zero.
    // -------------------------------------------------------------------------

    class FakeHalfYArchive final : public EffectiveSourceArchive {
    public:
        SourceMetadata metadata() const override {
            SourceMetadata meta;
            meta.rp = Real(10);
            meta.B  = Real(2);
            meta.sym.negative_m_is_conjugate = true;
            meta.sym.y_storage = YStorage::HalfWithEvenReflection;
            return meta;
        }

    protected:
        ModeData load_nonnegative_mode_patch_raw_(int m, PatchSide patch) const override {
            ModeData out;
            out.m = m;
            out.patch = patch;

            if (patch == PatchSide::Left) {
                out.X.values = {Real(-2), Real(-1)};
            } else {
                out.X.values = {Real(1), Real(2)};
            }

            out.Y.values = {Real(0), Real(1), Real(2)};

            const size_t nx = out.X.size();
            const size_t ny = out.Y.size();
            const size_t N  = nx * ny;

            for (auto& comp : out.data) {
                comp.assign(N, Complex(Real(0), Real(0)));
            }

            auto& f = out.component(Component::ll);

            for (size_t i = 0; i < nx; ++i) {
                const Real X = out.X.values[i];
                for (size_t j = 0; j < ny; ++j) {
                    const Real Y = out.Y.values[j];
                    const Real re = Real(100 * m) + Real(10) * X + Y;
                    const Real im = Real(1000 * m) + Real(100) * X + Real(10) * Y;
                    f[i * ny + j] = Complex(re, im);
                }
            }

            return out;
        }
    };

    Complex expected_ll_full(int m, Real X, Real Yabs) {
        return Complex(
                Real(100 * m) + Real(10) * X + Yabs,
                Real(1000 * m) + Real(100) * X + Real(10) * Yabs
        );
    }

    void test_half_y_expansion_right_patch() {
        FakeHalfYArchive archive;
        ModeData mode = archive.load_mode_patch(2, PatchSide::Right);

        expect_true(mode.patch == PatchSide::Right,
                    "test_half_y_expansion_right_patch: wrong patch");
        expect_true(mode.nx() == 2, "test_half_y_expansion_right_patch: nx should be 2");
        expect_true(mode.ny() == 5, "test_half_y_expansion_right_patch: ny should be 5");

        const std::vector<Real> expectedY = {
                Real(-2), Real(-1), Real(0), Real(1), Real(2)
        };

        for (size_t j = 0; j < expectedY.size(); ++j) {
            expect_real_eq(mode.Y.values[j], expectedY[j],
                           "test_half_y_expansion_right_patch: Y axis mismatch");
        }

        const auto& f = mode.component(Component::ll);

        // X = 1 row
        expect_complex_eq(f[0 * 5 + 0], expected_ll_full(2, Real(1), Real(2)),
                          "test_half_y_expansion_right_patch: f(1,-2)");
        expect_complex_eq(f[0 * 5 + 2], expected_ll_full(2, Real(1), Real(0)),
                          "test_half_y_expansion_right_patch: f(1,0)");
        expect_complex_eq(f[0 * 5 + 4], expected_ll_full(2, Real(1), Real(2)),
                          "test_half_y_expansion_right_patch: f(1,2)");

        // X = 2 row
        expect_complex_eq(f[1 * 5 + 0], expected_ll_full(2, Real(2), Real(2)),
                          "test_half_y_expansion_right_patch: f(2,-2)");
        expect_complex_eq(f[1 * 5 + 2], expected_ll_full(2, Real(2), Real(0)),
                          "test_half_y_expansion_right_patch: f(2,0)");
        expect_complex_eq(f[1 * 5 + 4], expected_ll_full(2, Real(2), Real(2)),
                          "test_half_y_expansion_right_patch: f(2,2)");
    }

    void test_negative_m_conjugation_right_patch() {
        FakeHalfYArchive archive;

        ModeData pos = archive.load_mode_patch(2, PatchSide::Right);
        ModeData neg = archive.load_mode_patch(-2, PatchSide::Right);

        const auto& fp = pos.component(Component::ll);
        const auto& fn = neg.component(Component::ll);

        expect_true(fp.size() == fn.size(),
                    "test_negative_m_conjugation_right_patch: size mismatch");

        for (size_t k = 0; k < fp.size(); ++k) {
            const Complex expected(fp[k].real(), -fp[k].imag());
            expect_complex_eq(fn[k], expected,
                              "test_negative_m_conjugation_right_patch: conjugation mismatch");
        }
    }

    void test_sample_native_right_patch() {
        FakeHalfYArchive archive;
        SourceMetadata meta = archive.metadata();
        ModeData mode = archive.load_mode_patch(2, PatchSide::Right);

        EffectiveSourceSampler sampler(meta, mode);

        // Interpolate in cell X in [1,2], Y in [0,1]
        const Complex got1 = sampler.sample_native(Component::ll, Real(1.5), Real(0.5));
        const Complex exp1 = expected_ll_full(2, Real(1.5), Real(0.5));
        expect_complex_eq(got1, exp1, "test_sample_native_right_patch: sample at (1.5,0.5)");

        const Complex got2 = sampler.sample_native(Component::ll, Real(1.5), Real(-0.5));
        const Complex exp2 = expected_ll_full(2, Real(1.5), Real(0.5));
        expect_complex_eq(got2, exp2, "test_sample_native_right_patch: sample at (1.5,-0.5)");
    }

    void test_sample_physical_right_patch() {
        FakeHalfYArchive archive;
        SourceMetadata meta = archive.metadata();
        ModeData mode = archive.load_mode_patch(2, PatchSide::Right);

        EffectiveSourceSampler sampler(meta, mode);

        // Want X = 1.5, Y = 0.5
        // X = (r - rp)/B, with rp=10, B=2 => r = 13
        // Y = -rp z => z = -0.05
        const Real r = Real(13);
        const Real z = Real(-0.05);

        const Complex got = sampler.sample_physical(Component::ll, r, z);
        const Complex exp = expected_ll_full(2, Real(1.5), Real(0.5));

        expect_complex_eq(got, exp, "test_sample_physical_right_patch");
    }

    void test_clamping_right_patch() {
        FakeHalfYArchive archive;
        SourceMetadata meta = archive.metadata();
        ModeData mode = archive.load_mode_patch(2, PatchSide::Right);

        EffectiveSourceSampler sampler(meta, mode);

        // Clamp x -> 1, y -> 2
        const Complex got = sampler.sample_native(Component::ll, Real(-100), Real(999));
        const Complex exp = expected_ll_full(2, Real(1), Real(2));

        expect_complex_eq(got, exp, "test_clamping_right_patch");
    }

    void test_left_patch_validation_and_sampling() {
        FakeHalfYArchive archive;
        SourceMetadata meta = archive.metadata();
        ModeData mode = archive.load_mode_patch(2, PatchSide::Left);

        expect_true(mode.patch == PatchSide::Left,
                    "test_left_patch_validation_and_sampling: wrong patch");
        expect_true(mode.X.max() < Real(0),
                    "test_left_patch_validation_and_sampling: left patch should have X<0");

        EffectiveSourceSampler sampler(meta, mode);

        const Complex got = sampler.sample_native(Component::ll, Real(-1.5), Real(0.5));
        const Complex exp = expected_ll_full(2, Real(-1.5), Real(0.5));

        expect_complex_eq(got, exp, "test_left_patch_validation_and_sampling");
    }

} // namespace

int main() {
    try {
        test_half_y_expansion_right_patch();
        test_negative_m_conjugation_right_patch();
        test_sample_native_right_patch();
        test_sample_physical_right_patch();
        test_clamping_right_patch();
        test_left_patch_validation_and_sampling();

        std::cout << "All EffectiveSource tests passed.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "EffectiveSource test FAILED: " << e.what() << '\n';
        return 1;
    }
}