//
// Created by Peter Zimmerman on 26.10.25.
//

#ifndef GHZ_NUMERIC_VECTORSGHZ_HPP
#define GHZ_NUMERIC_VECTORSGHZ_HPP

#pragma once
#ifndef GHZ_NUMERIC_VECTORS_HPP
#define GHZ_NUMERIC_VECTORS_HPP

#include "TeukTypes.hpp"
#include <array>
#include <cassert>

namespace ghz {
    using teuk::Real;
    using teuk::Complex;

/**
 * @brief 4-vector with scalar components (Real precision).
 */
    struct Vector4 {
        std::array<Real, 4> data{};

        constexpr Vector4(Real V0 = 0, Real V1 = 0, Real V2 = 0, Real V3 = 0) noexcept
                : data{V0, V1, V2, V3} {}

        constexpr Real& operator[](std::size_t i) noexcept {
            assert(i < 4);
            return data[i];
        }
        constexpr Real operator[](std::size_t i) const noexcept {
            assert(i < 4);
            return data[i];
        }

        [[nodiscard]] constexpr Vector4 operator+(const Vector4& o) const noexcept {
            return {data[0] + o.data[0], data[1] + o.data[1],
                    data[2] + o.data[2], data[3] + o.data[3]};
        }

        [[nodiscard]] constexpr Vector4 operator-(const Vector4& o) const noexcept {
            return {data[0] - o.data[0], data[1] - o.data[1],
                    data[2] - o.data[2], data[3] - o.data[3]};
        }

        [[nodiscard]] constexpr Vector4 operator*(Real s) const noexcept {
            return {s * data[0], s * data[1], s * data[2], s * data[3]};
        }

        constexpr Vector4& operator+=(const Vector4& o) noexcept {
            for (int i = 0; i < 4; ++i) data[i] += o.data[i];
            return *this;
        }

        constexpr Vector4& operator-=(const Vector4& o) noexcept {
            for (int i = 0; i < 4; ++i) data[i] -= o.data[i];
            return *this;
        }

        constexpr Vector4& operator*=(Real s) noexcept {
            for (auto& x : data) x *= s;
            return *this;
        }
    };

/**
 * @brief Complex 4-vector.
 */
    struct CVector4 {
        std::array<Complex, 4> data{};

        constexpr CVector4(Complex V0 = 0, Complex V1 = 0,
                           Complex V2 = 0, Complex V3 = 0) noexcept
                : data{V0, V1, V2, V3} {}

        constexpr Complex& operator[](std::size_t i) noexcept {
            assert(i < 4);
            return data[i];
        }
        constexpr Complex operator[](std::size_t i) const noexcept {
            assert(i < 4);
            return data[i];
        }

        [[nodiscard]] constexpr CVector4 operator*(Complex s) const noexcept {
            return {s * data[0], s * data[1], s * data[2], s * data[3]};
        }

        [[nodiscard]] CVector4 conj() const noexcept {
            return {std::conj(data[0]), std::conj(data[1]),
                    std::conj(data[2]), std::conj(data[3])};
        }
    };


/**
 * @brief Real symmetric 4×4 matrix (e.g. metric tensor).
 * Stores only 10 independent components.
 */
    struct SymmetricMatrix4 {
        std::array<Real, 10> data{}; // Storage order: (00, 01, 02, 03, 11, 12, 13, 22, 23, 33)

        constexpr SymmetricMatrix4() noexcept = default;

        constexpr SymmetricMatrix4(std::initializer_list<Real> vals) noexcept {
            assert(vals.size() == 10);
            std::copy(vals.begin(), vals.end(), data.begin());
        }

        // Access symmetric components
        constexpr Real& operator()(int i, int j) noexcept {
            assert(i >= 0 && i < 4 && j >= 0 && j < 4);
            return data[index(i, j)];
        }

        constexpr Real operator()(int i, int j) const noexcept {
            assert(i >= 0 && i < 4 && j >= 0 && j < 4);
            return data[index(i, j)];
        }

        [[nodiscard]] constexpr SymmetricMatrix4 operator+(const SymmetricMatrix4& o) const noexcept {
            SymmetricMatrix4 r;
            for (int k = 0; k < 10; ++k) r.data[k] = data[k] + o.data[k];
            return r;
        }

        [[nodiscard]] constexpr SymmetricMatrix4 operator-(const SymmetricMatrix4& o) const noexcept {
            SymmetricMatrix4 r;
            for (int k = 0; k < 10; ++k) r.data[k] = data[k] - o.data[k];
            return r;
        }

        [[nodiscard]] constexpr SymmetricMatrix4 operator*(Real s) const noexcept {
            SymmetricMatrix4 r;
            for (int k = 0; k < 10; ++k) r.data[k] = s * data[k];
            return r;
        }

        // Apply matrix to vector: v'ᵢ = gᵢⱼ vʲ
        [[nodiscard]] constexpr Vector4 operator*(const Vector4& v) const noexcept {
            Vector4 r;
            for (int i = 0; i < 4; ++i) {
                Real sum = 0;
                for (int j = 0; j < 4; ++j)
                    sum += (*this)(i, j) * v[j];
                r[i] = sum;
            }
            return r;
        }

        // Optional: return full 4×4 matrix view
        [[nodiscard]] constexpr std::array<std::array<Real, 4>, 4> full() const noexcept {
            std::array<std::array<Real, 4>, 4> m{};
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    m[i][j] = (*this)(i, j);
            return m;
        }

    private:
        // Map (i,j) with i ≤ j to compact index
        static constexpr int index(int i, int j) noexcept {
            if (i > j) std::swap(i, j);
            constexpr int map[4][4] = {
                    {0, 1, 2, 3},
                    {1, 4, 5, 6},
                    {2, 5, 7, 8},
                    {3, 6, 8, 9}
            };
            return map[i][j];
        }
    };

} // namespace ghz

#endif // GHZ_NUMERIC_VECTORS_HPP


#endif //GHZ_NUMERIC_VECTORSGHZ_HPP
