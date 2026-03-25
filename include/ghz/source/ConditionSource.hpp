//
// Created by Peter Zimmerman on 17.03.26.
//
#ifndef GHZ_SOURCE_CONDITIONSOURCE_HPP
#define GHZ_SOURCE_CONDITIONSOURCE_HPP
#pragma once

#include "ghz/source/EffectiveSource.hpp"

#include <array>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace ghz::source {

    class ConditionSource {
    public:
        using Real = teuk::Real;
        using Complex = teuk::Complex;

        ConditionSource(const EffectiveSourceArchive& archive,
                        int m,
                        PatchSide patch);

        // Precompute source values + cubic slopes on a fixed solver z-slice.
        void set_z_slice(Real z_fixed);

        [[nodiscard]] bool has_active_z_slice() const noexcept { return has_z_slice_; }
        [[nodiscard]] Real z_slice() const {
            if (!has_z_slice_) {
                throw std::runtime_error("ConditionSource::z_slice(): z slice not set");
            }
            return z_fixed_;
        }

        [[nodiscard]] int m() const noexcept { return m_; }
        [[nodiscard]] PatchSide patch() const noexcept { return patch_; }

        [[nodiscard]] Real r_min() const noexcept;
        [[nodiscard]] Real r_max() const noexcept;
        [[nodiscard]] bool contains_r(Real r) const noexcept;

        [[nodiscard]] Complex eval(Component c, Real r) const;
        [[nodiscard]] Complex operator()(Component c, Real r) const {
            return eval(c, r);
        }

        [[nodiscard]] std::vector<Complex>
        eval_on_r_grid(Component c, const std::vector<Real>& r_nodes) const;

        [[nodiscard]] Complex Tll(Real r) const { return eval(Component::ll, r); }
        [[nodiscard]] Complex Tlm(Real r) const { return eval(Component::lm, r); }
        [[nodiscard]] Complex Tln(Real r) const { return eval(Component::ln, r); }

        [[nodiscard]] std::vector<Complex> Tll_on_r_grid(const std::vector<Real>& r_nodes) const {
            return eval_on_r_grid(Component::ll, r_nodes);
        }

        [[nodiscard]] std::vector<Complex> Tlm_on_r_grid(const std::vector<Real>& r_nodes) const {
            return eval_on_r_grid(Component::lm, r_nodes);
        }

        [[nodiscard]] std::vector<Complex> Tln_on_r_grid(const std::vector<Real>& r_nodes) const {
            return eval_on_r_grid(Component::ln, r_nodes);
        }

    private:
        SourceMetadata meta_;
        EffectiveSourceSampler sampler_;
        int m_;
        PatchSide patch_;

        // Native X nodes of this patch
        std::vector<Real> x_nodes_;

        // Precomputed values and cubic slopes on current fixed z-slice
        std::array<std::vector<Complex>, component_count> slice_values_{};
        std::array<std::vector<Complex>, component_count> slice_slopes_{};

        bool has_z_slice_ = false;
        Real z_fixed_ = Real(0);

        static size_t lower_cell_(const std::vector<Real>& grid, Real x);

        [[nodiscard]] Complex cubic_interp_(const std::vector<Complex>& f,
                                            const std::vector<Complex>& d,
                                            Real Xq) const;

        static std::vector<Real> pchip_slopes_real_(const std::vector<Real>& x,
                                                    const std::vector<Real>& y);

        static std::vector<Complex> pchip_slopes_complex_(const std::vector<Real>& x,
                                                          const std::vector<Complex>& y);

        static bool same_sign_nonzero_(Real a, Real b) noexcept;
        static Real abs_(Real x) noexcept;
    };

} // namespace ghz::source

#endif // GHZ_SOURCE_CONDITIONSOURCE_HPP