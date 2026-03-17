//
// Created by Peter Zimmerman on 04.03.26.
//

#ifndef GHZ_SOURCE_CONDITIONSOURCE_HPP
#define GHZ_SOURCE_CONDITIONSOURCE_HPP
#pragma once

#include "ghz/core/GhzTypes.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cmath>

namespace ghz::source {

    using teuk::Real;
    using teuk::Complex;


    class ConditionSource {
    public:
        using Real = teuk::Real;
        using Complex = teuk::Complex;

        struct UniformGrid {
            size_t Nr, Nz;
            Real r_min, r_max;
            Real z_min, z_max;
        };

        enum class FileFormat {
            BinaryComplexInterleaved,
            TextComplexTwoCols,
            TextRealOnly // optional; imag=0
        };

        struct WindowParams {
            Real center_r;
            Real halfwidth_r;
        };

        ConditionSource(const UniformGrid& grid);

        void set_window_params(const WindowParams& w);

        void load_from_file(const std::string& path, FileFormat fmt);

        void apply_window_inplace();

        std::vector<Complex>
        interpolate_to_collocation(const std::vector<Real>& r_nodes,
                                   const std::vector<Real>& z_nodes) const;

    private:
        UniformGrid grid_;
        WindowParams window_;

        std::vector<Complex> data_;

        Real window_value_(Real r, Real z) const;
        Complex bilinear_(Real r, Real z) const;
    };

}
#endif //GHZ_SOURCE_CONDITIONSOURCE_HPP
