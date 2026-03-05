//
// Created by Peter Zimmerman on 04.03.26.
//
#pragma once
/**
 * @brief Read vectorized effective-source data (uniform r,z grid, row-major) and apply a smooth window,
 *        then interpolate onto spectral collocation grids (Chebyshev in r, LGL in z). \n
 *
 * Assumptions / conventions:\n
 * - Source data lives on a *uniform* tensor grid (r_i, z_j), with sizes (Nr_u, Nz_u). \n
 * - Stored in row-major order: idx(i,j) = i*Nz_u + j, where i is r-index, j is z-index. \n
 * - Assumes one  file per mode (m,kr,kz)  \n
 * - Interpolation provided here: bilinear (robust and fast). You can swap in higher-order later. \n
 *
 */

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <type_traits>
#include <optional>
#include "ghz/source/ConditionSource.hpp"



namespace ghz {

    ConditionSource::ConditionSource(const UniformGrid& grid)
            : grid_(grid)
    {
        data_.resize(grid_.Nr * grid_.Nz);
    }

    void ConditionSource::set_window_params(const WindowParams& w)
    {
        window_ = w;
    }

    void ConditionSource::load_from_file(const std::string& path, FileFormat fmt)
    {
        const size_t N = grid_.Nr * grid_.Nz;

        std::ifstream in;

        if (fmt == FileFormat::BinaryComplexInterleaved) {

            in.open(path, std::ios::binary);
            if (!in)
                throw std::runtime_error("ConditionSource: cannot open file " + path);

            std::vector<Real> buffer(2 * N);

            in.read(reinterpret_cast<char*>(buffer.data()),
                    static_cast<std::streamsize>(buffer.size() * sizeof(Real)));

            if (!in)
                throw std::runtime_error("ConditionSource: binary read failed");

            for (size_t k = 0; k < N; ++k)
                data_[k] = Complex(buffer[2*k], buffer[2*k + 1]);
        }

        else if (fmt == FileFormat::TextComplexTwoCols) {

            in.open(path);
            if (!in)
                throw std::runtime_error("ConditionSource: cannot open file " + path);

            Real re, im;

            for (size_t k = 0; k < N; ++k) {

                if (!(in >> re >> im))
                    throw std::runtime_error("ConditionSource: unexpected EOF");

                data_[k] = Complex(re, im);
            }
        }

        else if (fmt == FileFormat::TextRealOnly) {

            in.open(path);
            if (!in)
                throw std::runtime_error("ConditionSource: cannot open file " + path);

            Real re;

            for (size_t k = 0; k < N; ++k) {

                if (!(in >> re))
                    throw std::runtime_error("ConditionSource: unexpected EOF");

                data_[k] = Complex(re, 0.0);
            }
        }
    }

    void ConditionSource::apply_window_inplace()
    {
        const Real dr = (grid_.r_max - grid_.r_min) / Real(grid_.Nr - 1);
        const Real dz = (grid_.z_max - grid_.z_min) / Real(grid_.Nz - 1);

        for (size_t i = 0; i < grid_.Nr; ++i) {

            Real r = grid_.r_min + Real(i) * dr;

            for (size_t j = 0; j < grid_.Nz; ++j) {

                Real z = grid_.z_min + Real(j) * dz;

                Real W = window_value_(r, z);

                data_[i * grid_.Nz + j] *= W;
            }
        }
    }

    ConditionSource::Real
    ConditionSource::window_value_(Real r, Real z) const
    {
        const Real c = window_.center_r;
        const Real w = window_.halfwidth_r;

        if (w <= 0.0)
            return 1.0;

        Real s = (r - c) / w;

        if (std::abs(s) >= 1.0)
            return 0.0;

        Real t = 1.0 - s*s;

        return std::exp(-1.0 / t);  // C^∞ bump
    }

    ConditionSource::Complex
    ConditionSource::bilinear_(Real r, Real z) const
    {
        const Real dr = (grid_.r_max - grid_.r_min) / Real(grid_.Nr - 1);
        const Real dz = (grid_.z_max - grid_.z_min) / Real(grid_.Nz - 1);

        r = std::clamp(r, grid_.r_min, grid_.r_max);
        z = std::clamp(z, grid_.z_min, grid_.z_max);

        Real fr = (r - grid_.r_min) / dr;
        Real fz = (z - grid_.z_min) / dz;

        size_t i0 = static_cast<size_t>(std::floor(fr));
        size_t j0 = static_cast<size_t>(std::floor(fz));

        if (i0 >= grid_.Nr - 1) i0 = grid_.Nr - 2;
        if (j0 >= grid_.Nz - 1) j0 = grid_.Nz - 2;

        size_t i1 = i0 + 1;
        size_t j1 = j0 + 1;

        Real tr = fr - Real(i0);
        Real tz = fz - Real(j0);

        const Complex& f00 = data_[i0 * grid_.Nz + j0];
        const Complex& f10 = data_[i1 * grid_.Nz + j0];
        const Complex& f01 = data_[i0 * grid_.Nz + j1];
        const Complex& f11 = data_[i1 * grid_.Nz + j1];

        Complex a = (1.0 - tr) * f00 + tr * f10;
        Complex b = (1.0 - tr) * f01 + tr * f11;

        return (1.0 - tz) * a + tz * b;
    }

    std::vector<ConditionSource::Complex>
    ConditionSource::interpolate_to_collocation(
            const std::vector<Real>& r_nodes,
            const std::vector<Real>& z_nodes) const
    {
        const size_t Nr = r_nodes.size();
        const size_t Nz = z_nodes.size();

        std::vector<Complex> out(Nr * Nz);

        for (size_t i = 0; i < Nr; ++i) {

            for (size_t j = 0; j < Nz; ++j) {

                out[i * Nz + j] = bilinear_(r_nodes[i], z_nodes[j]);
            }
        }

        return out;
    }

} // namespace ghz