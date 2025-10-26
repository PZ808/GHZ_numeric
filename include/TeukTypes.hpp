//
// Created by Peter Zimmerman on 26.10.25.
//

#ifndef GHZ_NUMERIC_TEUKTYPES_HPP
#define GHZ_NUMERIC_TEUKTYPES_HPP

#pragma once
/**
 * @file TeukTypes.hpp
 * @brief Type and precision definitions for the GHZ_numeric project.
 */

#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>
#include <numbers>

namespace teuk {

/**
 * @brief Enumeration for selecting floating-point precision at compile time.
 */
    enum class PrecisionLevel : int {
        Float       = 0,
        Double      = 1,
        LongDouble  = 2
    };

/**
 * @brief Choose precision level (can be overridden by compiler definition).
 * Example: -TEUK_PRECISION=2 for long double
 */
#ifndef TEUK_PRECISION
#  define TEUK_PRECISION 1
#endif

// Primary template — not defined.
    template <PrecisionLevel> struct PrecisionTraits;

// Specializations:
    template <> struct PrecisionTraits<PrecisionLevel::Float> {
        using Real       = float;
        using LongReal   = double;
        static constexpr auto epsilon = std::numeric_limits<Real>::epsilon();
        static constexpr int mantissa_digits = std::numeric_limits<Real>::digits;
    };

    template <> struct PrecisionTraits<PrecisionLevel::Double> {
        using Real       = double;
        using LongReal   = long double;
        static constexpr auto epsilon = std::numeric_limits<Real>::epsilon();
        static constexpr int mantissa_digits = std::numeric_limits<Real>::digits;
    };

    template <> struct PrecisionTraits<PrecisionLevel::LongDouble> {
        using Real       = long double;
        using LongReal   = long double;
        static constexpr auto epsilon = std::numeric_limits<Real>::epsilon();
        static constexpr int mantissa_digits = std::numeric_limits<Real>::digits;
    };

// Alias for the active precision
    using ActivePrecision = PrecisionTraits<
            static_cast<PrecisionLevel>(TEUK_PRECISION)
    >;

    using Real       = typename ActivePrecision::Real;
    using LongReal   = typename ActivePrecision::LongReal;
    using Complex    = std::complex<Real>;
    using LongComplex= std::complex<LongReal>;

// Imaginary unit
    inline constexpr Complex I{0, 1};

// Complex utilities
    inline constexpr Complex CNAN{std::numeric_limits<Real>::quiet_NaN(),
                                  std::numeric_limits<Real>::quiet_NaN()};

    template <typename T>
    inline constexpr bool Cisnan(const std::complex<T>& z) noexcept {
        return std::isnan(std::real(z)) || std::isnan(std::imag(z));
    }

    template <typename T>
    inline constexpr T Cmaxpart(const std::complex<T>& z) noexcept {
        return std::fmax(std::fabs(std::real(z)), std::fabs(std::imag(z)));
    }

// Precision-dependent math wrappers (constexpr)
    template <typename T>
    inline constexpr T Fmax(T a, T b) noexcept { return std::fmax(a, b); }

    template <typename T>
    inline constexpr T Fmin(T a, T b) noexcept { return std::fmin(a, b); }

    template <typename T>
    inline constexpr T Fabs(T x) noexcept { return std::fabs(x); }

    template <typename T>
    inline constexpr T Sqrt(T x) noexcept { return std::sqrt(x); }

} // namespace teuk


#endif //GHZ_NUMERIC_TEUKTYPES_HPP
