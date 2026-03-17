//
// Created by Peter Zimmerman on 17.03.26.
//

#ifndef GHZ_SOURCE_EFFECTIVESOURCE_HPP
#define GHZ_SOURCE_EFFECTIVESOURCE_HPP
#pragma once

#include "ghz/core/GhzTypes.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>



namespace ghz::source {

    using teuk::Real;
    using teuk::Complex;

    enum class PatchSide : std::uint32_t {
        Left  = 0,
        Right = 1
    };
    inline const char* patch_tag(PatchSide p) noexcept {
        return p == PatchSide::Left ? "left" : "right";
    }

    // -------------------------------------------------------------------------
    // Axis
    // -------------------------------------------------------------------------

    struct Axis {
        std::vector<Real> values;

        size_t size() const noexcept { return values.size(); }
        bool empty() const noexcept { return values.empty(); }

        Real min() const {
            if (values.empty()) {
                throw std::runtime_error("Axis::min(): empty axis");
            }
            return values.front();
        }

        Real max() const {
            if (values.empty()) {
                throw std::runtime_error("Axis::max(): empty axis");
            }
            return values.back();
        }

        bool is_uniform(Real tol = Real(1e-12)) const {
            if (values.size() < 3) return true;

            const Real h = values[1] - values[0];
            for (size_t i = 1; i + 1 < values.size(); ++i) {
                if (teuk::Fabs((values[i + 1] - values[i]) - h) > tol) {
                    return false;
                }
            }
            return true;
        }

        Real dx_uniform() const {
            if (values.size() < 2) {
                throw std::runtime_error("Axis::dx_uniform(): axis too short");
            }
            return (values.back() - values.front()) / Real(values.size() - 1);
        }
    };

    // -------------------------------------------------------------------------
    // Components
    // -------------------------------------------------------------------------

    enum class Component : std::uint32_t {
        ll = 0,
        ln,
        lm,
        lmb,
        nn,
        nm,
        nmb,
        mm,
        mmb,
        count
    };

    constexpr size_t component_index(Component c) noexcept {
        return static_cast<size_t>(c);
    }

    constexpr size_t component_count = component_index(Component::count);

    // -------------------------------------------------------------------------
    // Storage symmetry / metadata
    // -------------------------------------------------------------------------

    enum class YStorage : std::uint32_t {
        Full = 0,
        HalfWithEvenReflection = 1
    };
    // -------------------------------------------------------------------------
    // Symmetry / metadata
    // -------------------------------------------------------------------------
    struct SourceSymmetry {
        // If true, raw files only need to store m >= 0 and negative m
        // is reconstructed as complex conjugation of +|m|.
        bool negative_m_is_conjugate = true;

        // How the Y axis is stored on disk.
        // Full                   : Y contains the entire domain.
        // HalfWithEvenReflection : raw data stores Y >= 0 including Y = 0,
        //                          and load_mode_patch() expands to full Y.
        YStorage y_storage = YStorage::Full;
    };

    struct SourceMetadata {
        Real rp{};
        Real B{};
        SourceSymmetry sym{};
    };

    // -------------------------------------------------------------------------
    // Mode data for a single m
    //
    // Flattening convention:
    //   data[c][i * ny + j]  <=>  component c at (X_i, Y_j)
    // -------------------------------------------------------------------------

    struct ModeData {
        int m = 0;
        PatchSide patch = PatchSide::Left;
        Axis X;
        Axis Y;

        std::array<std::vector<Complex>, component_count> data;

        size_t nx() const noexcept { return X.size(); }
        size_t ny() const noexcept { return Y.size(); }

        size_t flat_index(size_t i, size_t j) const noexcept {
            return i * ny() + j;
        }

        const std::vector<Complex>& component(Component c) const {
            return data[component_index(c)];
        }

        std::vector<Complex>& component(Component c) {
            return data[component_index(c)];
        }

        void validate() const {
            validate_axis_strictly_increasing(X, "ModeData::X");
            validate_axis_strictly_increasing(Y, "ModeData::Y");

            const size_t N = nx() * ny();
            for (size_t k = 0; k < data.size(); ++k) {
                if (data[k].size() != N) {
                    throw std::runtime_error("ModeData::validate(): component has wrong size");
                }
            }
        }

        static void validate_axis_strictly_increasing(const Axis& a,
                                                      const std::string& name) {
            if (a.empty()) {
                throw std::runtime_error(name + ": empty axis");
            }
            for (size_t i = 1; i < a.size(); ++i) {
                if (!(a.values[i] > a.values[i - 1])) {
                    throw std::runtime_error(name + ": axis must be strictly increasing");
                }
            }
        }

        static void validate_raw_against_metadata(const ModeData& mode,
                                                  const SourceMetadata& meta) {
            mode.validate();

            switch (meta.sym.y_storage) {
                case YStorage::Full:
                    break;

                case YStorage::HalfWithEvenReflection:
                    if (mode.Y.empty()) {
                        throw std::runtime_error(
                                "ModeData::validate_raw_against_metadata(): half-Y storage with empty Y axis");
                    }
                    if (mode.Y.min() < Real(0)) {
                        throw std::runtime_error(
                                "ModeData::validate_raw_against_metadata(): half-Y storage requires Y >= 0");
                    }
                    if (mode.Y.values.front() != Real(0)) {
                        throw std::runtime_error(
                                "ModeData::validate_raw_against_metadata(): half-Y storage must include Y = 0 as first point");
                    }
                    break;

                default:
                    throw std::runtime_error(
                            "ModeData::validate_raw_against_metadata(): unsupported YStorage");
            }

            if (mode.patch == PatchSide::Left) {
                if (!(mode.X.max() < Real(0))) {
                    throw std::runtime_error(
                            "ModeData::validate_raw_against_metadata(): left patch must satisfy X < 0");
                }
            } else {
                if (!(mode.X.min() > Real(0))) {
                    throw std::runtime_error(
                            "ModeData::validate_raw_against_metadata(): right patch must satisfy X > 0");
                }
            }
        }
    };

    // -------------------------------------------------------------------------
    // Archive interface
    //
    // Contract for derived classes:
    //   - metadata() returns archive-wide metadata.
    //   - load_nonnegative_mode_patch_raw_(m, patch) loads the stored raw mode
    //     for m >= 0 on the specified radial patch.
    // ...
    // Public load_mode_patch(m, patch):
    //   - raw data may be stored with full Y or half Y depending on metadata().
    //   - handles m < 0 via conjugation if allowed
    //   - expands half-Y storage to a canonical full-Y ModeData
    //   - validates the result
    // -------------------------------------------------------------------------

    class EffectiveSourceArchive {
    public:
        virtual ~EffectiveSourceArchive() = default;

        virtual SourceMetadata metadata() const = 0;

        ModeData load_mode_patch(int m, PatchSide patch) const {
            const SourceMetadata meta = metadata();

            if (m >= 0) {
                ModeData raw = load_nonnegative_mode_patch_raw_(m, patch);
                raw.m = m;
                raw.patch = patch;

                ModeData::validate_raw_against_metadata(raw, meta);

                ModeData expanded = expand_y_if_needed_(meta, std::move(raw));
                expanded.m = m;
                expanded.patch = patch;
                expanded.validate();
                return expanded;
            }

            if (!meta.sym.negative_m_is_conjugate) {
                throw std::runtime_error(
                        "EffectiveSourceArchive::load_mode_patch(): negative m requested but symmetry is unavailable");
            }

            ModeData raw = load_nonnegative_mode_patch_raw_(-m, patch);
            raw.m = m;
            raw.patch = patch;

            ModeData::validate_raw_against_metadata(raw, meta);

            for (auto& comp : raw.data) {
                for (auto& z : comp) {
                    z = Complex(z.real(), -z.imag());
                }
            }

            ModeData expanded = expand_y_if_needed_(meta, std::move(raw));
            expanded.m = m;
            expanded.patch = patch;
            expanded.validate();
            return expanded;
        }

    protected:
        virtual ModeData load_nonnegative_mode_patch_raw_(int m, PatchSide patch) const = 0;

    private:
        static ModeData expand_y_if_needed_(const SourceMetadata& meta,
                                            ModeData mode) {
            if (meta.sym.y_storage == YStorage::Full) {
                return mode;
            }

            if (meta.sym.y_storage != YStorage::HalfWithEvenReflection) {
                throw std::runtime_error(
                        "EffectiveSourceArchive::expand_y_if_needed_(): "
                        "unsupported YStorage"
                );
            }

            const size_t nx  = mode.nx();
            const size_t nyh = mode.ny();

            if (nyh == 0) {
                throw std::runtime_error(
                        "EffectiveSourceArchive::expand_y_if_needed_(): empty half-Y axis"
                );
            }

            ModeData out;
            out.m = mode.m;
            out.X = mode.X;

            out.patch = mode.patch;
            out.Y.values.clear();
            out.Y.values.reserve(2 * nyh - 1);

            // negative side in increasing order: [-y_last, ..., -y_2, -y_1]
            for (size_t j = nyh - 1; j >= 1; --j) {
                out.Y.values.push_back(-mode.Y.values[j]);
                if (j == 1) break; // avoid size_t underflow
            }

            // nonnegative side: [0, y1, y2, ...]
            for (const Real y : mode.Y.values) {
                out.Y.values.push_back(y);
            }

            const size_t nyf = out.Y.size();

            for (size_t c = 0; c < component_count; ++c) {
                out.data[c].resize(nx * nyf);

                for (size_t i = 0; i < nx; ++i) {

                    // left half: negative Y
                    for (size_t j = 1; j < nyh; ++j) {
                        const size_t j_left = nyh - 1 - j;
                        out.data[c][i * nyf + j_left] =
                                mode.data[c][i * nyh + j];
                    }

                    // right half: 0 and positive Y
                    for (size_t j = 0; j < nyh; ++j) {
                        const size_t j_right = (nyh - 1) + j;
                        out.data[c][i * nyf + j_right] =
                                mode.data[c][i * nyh + j];
                    }
                }
            }

            return out;
        } // expand_y_if_needed_()
    };



    // -------------------------------------------------------------------------
    // EffectiveSourceSampler
    //
    // Assumes ModeData is already in canonical full-Y form, e.g. returned by
    // EffectiveSourceArchive::load_mode_patch().
    // -------------------------------------------------------------------------

    class EffectiveSourceSampler {
    public:
        EffectiveSourceSampler(SourceMetadata meta, ModeData data)
                : meta_(std::move(meta)), data_(std::move(data))
        {
            data_.validate();
        }

        Complex sample_native(Component c, Real Xq, Real Yq) const {
            return bilinear_(data_.component(c), Xq, Yq);
        }

        // Native-to-physical map:
        //   X = (r - rp) / B
        //   Y = -rp z
        Complex sample_physical(Component c, Real r, Real z) const {
            const Real Xq = (r - meta_.rp) / meta_.B;
            const Real Yq = -meta_.rp * z;
            return sample_native(c, Xq, Yq);
        }

        const SourceMetadata& metadata() const noexcept { return meta_; }
        const ModeData& mode_data() const noexcept { return data_; }

    private:
        SourceMetadata meta_;
        ModeData data_;

        Complex bilinear_(const std::vector<Complex>& f, Real x, Real y) const {
            const Axis& X = data_.X;
            const Axis& Y = data_.Y;

            if (X.size() < 2 || Y.size() < 2) {
                throw std::runtime_error(
                        "EffectiveSourceSampler::bilinear_(): axes too short"
                );
            }

            x = std::clamp(x, X.min(), X.max());
            y = std::clamp(y, Y.min(), Y.max());

            const size_t i0 = lower_cell_(X.values, x);
            const size_t j0 = lower_cell_(Y.values, y);
            const size_t i1 = i0 + 1;
            const size_t j1 = j0 + 1;

            const Real x0 = X.values[i0];
            const Real x1 = X.values[i1];
            const Real y0 = Y.values[j0];
            const Real y1 = Y.values[j1];

            const Real tx = (x1 != x0) ? (x - x0) / (x1 - x0) : Real(0);
            const Real ty = (y1 != y0) ? (y - y0) / (y1 - y0) : Real(0);

            const Complex& f00 = f[data_.flat_index(i0, j0)];
            const Complex& f10 = f[data_.flat_index(i1, j0)];
            const Complex& f01 = f[data_.flat_index(i0, j1)];
            const Complex& f11 = f[data_.flat_index(i1, j1)];

            const Complex a = (Real(1) - tx) * f00 + tx * f10;
            const Complex b = (Real(1) - tx) * f01 + tx * f11;

            return (Real(1) - ty) * a + ty * b;
        }

        static size_t lower_cell_(const std::vector<Real>& grid, Real x) {
            if (grid.size() < 2) {
                throw std::runtime_error(
                        "EffectiveSourceSampler::lower_cell_(): grid too short"
                );
            }

            auto it = std::upper_bound(grid.begin(), grid.end(), x);

            if (it == grid.begin()) {
                return 0;
            }

            if (it == grid.end()) {
                return grid.size() - 2;
            }

            return static_cast<size_t>(std::distance(grid.begin(), it) - 1);
        }
    };

    class TwoPatchEffectiveSource {
    public:
        TwoPatchEffectiveSource(SourceMetadata meta,
                                EffectiveSourceSampler left,
                                EffectiveSourceSampler right)
                : meta_(std::move(meta)),
                  left_(std::move(left)),
                  right_(std::move(right)) {}

        Complex sample_native(Component c, Real Xq, Real Yq) const {
            if (Xq < Real(0)) return left_.sample_native(c, Xq, Yq);
            if (Xq > Real(0)) return right_.sample_native(c, Xq, Yq);
            throw std::runtime_error("TwoPatchEffectiveSource::sample_native(): X=0 lies in excluded particle region");
        }

        Complex sample_physical(Component c, Real r, Real z) const {
            const Real Xq = (r - meta_.rp) / meta_.B;
            const Real Yq = -meta_.rp * z;
            return sample_native(c, Xq, Yq);
        }
        const EffectiveSourceSampler& left() const noexcept { return left_; }
        const EffectiveSourceSampler& right() const noexcept { return right_; }
    private:
        SourceMetadata meta_;
        EffectiveSourceSampler left_;
        EffectiveSourceSampler right_;
    };

} // namespace ghz::source

#endif // GHZ_SOURCE_EFFECTIVESOURCE_HPP