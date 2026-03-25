#ifndef GHZ_SOURCE_CONDITIONBOUNDARYDATA_HPP
#define GHZ_SOURCE_CONDITIONBOUNDARYDATA_HPP

#pragma once

#include "ghz/source/BinaryBoundaryDataArchive.hpp"

#include <array>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace ghz::source {

    class ConditionBoundaryData {
    public:
        using Real = teuk::Real;
        using Complex = teuk::Complex;

        ConditionBoundaryData(const BinaryBoundaryDataArchive& archive,
                              int m,
                              PatchSide patch);

        void set_z_slice(Real z_fixed);

        [[nodiscard]] bool has_active_z_slice() const noexcept { return has_z_slice_; }

        [[nodiscard]] Real z_slice() const {
            if (!has_z_slice_) {
                throw std::runtime_error("ConditionBoundaryData::z_slice(): z slice not set");
            }
            return z_fixed_;
        }

        [[nodiscard]] int m() const noexcept { return m_; }
        [[nodiscard]] PatchSide patch() const noexcept { return patch_; }

        [[nodiscard]] const BoundaryMetadata& metadata() const noexcept { return meta_; }

        [[nodiscard]] Real r_boundary() const noexcept { return meta_.r_boundary; }

        [[nodiscard]] Complex eval(BoundaryComponent c) const;
        [[nodiscard]] Complex operator()(BoundaryComponent c) const { return eval(c); }

        [[nodiscard]] Complex delta() const { return eval(BoundaryComponent::delta); }
        [[nodiscard]] Complex deltaPrime() const { return eval(BoundaryComponent::deltaPrime); }

    private:
        BoundaryMetadata meta_{};
        BoundaryModeData mode_data_{};
        int m_{0};
        PatchSide patch_{PatchSide::Left};

        std::vector<Real> y_nodes_;

        std::array<std::vector<Complex>, boundary_component_count> y_values_{};
        std::array<std::vector<Complex>, boundary_component_count> y_slopes_{};

        std::array<Complex, boundary_component_count> conditioned_values_{};

        bool has_z_slice_ = false;
        Real z_fixed_ = Real(0);

        static size_t lower_cell_(const std::vector<Real>& grid, Real x);

        [[nodiscard]] Complex cubic_interp_(const std::vector<Complex>& f,
                                            const std::vector<Complex>& d,
                                            Real yq) const;

        [[nodiscard]] Real canonical_y_(Real z) const;

        static std::vector<Real> pchip_slopes_real_(const std::vector<Real>& x,
                                                    const std::vector<Real>& y);

        static std::vector<Complex> pchip_slopes_complex_(const std::vector<Real>& x,
                                                          const std::vector<Complex>& y);

        static bool same_sign_nonzero_(Real a, Real b) noexcept;
        static Real abs_(Real x) noexcept;
    };

} // namespace ghz::source

#endif // GHZ_SOURCE_CONDITIONBOUNDARYDATA_HPP