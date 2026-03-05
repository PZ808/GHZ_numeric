//
// Created by Peter Zimmerman on 04.03.26.
//

#ifndef HS1DATA_DATA_CONDITIONSOURCE_HPP
#define HS1DATA_DATA_CONDITIONSOURCE_HPP

#pragma once
#include <vector>
#include <string>
#include <complex>

namespace ghz {

    class ConditionSource {
    public:
        using Real = double;
        using Complex = std::complex<Real>;

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
#endif //HS1DATA_DATA_CONDITIONSOURCE_HPP
