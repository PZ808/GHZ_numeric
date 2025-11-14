//
// Created by Peter Zimmerman on 25.10.25.
//
#ifndef GHZ_NUMERIC_MATHMACROS_HPP
#define GHZ_NUMERIC_MATHMACROS_HPP
#pragma once
/**
 * @file math_macros.hpp
 * @brief Safe, type-generic math utilities compatible with Boost.Multiprecision.
 *
 *  - Provides both C-style macros (for legacy code) and C++ inline constexpr functions.
 *  - All inline functions are precision-generic and work with teuk::Real / Complex types.
 */

#include <cmath>
#include <complex>
#include <type_traits>
#include <array>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "GhzTypes.hpp"  // for teuk::Real, Complex, etc.

// ======================================================
//  MACROS (legacy support — prefer inline functions!)
// ======================================================

#define SQR(x)    ((x) * (x))
#define CUBE(x)   ((x) * (x) * (x))
#define POW2(x)   SQR(x)
#define POW3(x)   CUBE(x)

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef CLAMP
#define CLAMP(v, lo, hi) (((v) < (lo)) ? (lo) : (((v) > (hi)) ? (hi) : (v)))
#endif
#ifndef SIGN
#define SIGN(x) (((x) > 0) - ((x) < 0))
#endif

// ======================================================
//  NAMESPACE math — safer, ADL-friendly alternatives
// ======================================================
namespace math {

// Helper for ADL to find boost::multiprecision overloads
    using std::abs;
    using std::sqrt;
    using std::sin;
    using std::cos;
    using std::tan;
    using std::exp;
    using std::log;
    using std::pow;

//----------------------------------------------------
// Simple arithmetic helpers
//----------------------------------------------------
    template <typename T>
    constexpr T sqr(const T& x) noexcept { return x * x; }

    template <typename T>
    constexpr T cube(const T& x) noexcept { return x * x * x; }

    template <typename T>
    constexpr T pow2(const T& x) noexcept { return x * x; }

    template <typename T>
    constexpr T pow3(const T& x) noexcept { return x * x * x; }

    template <typename T>
    constexpr T pow4(const T& x) noexcept { return (x * x) * (x * x); }

//----------------------------------------------------
// Min / Max / Clamp
//----------------------------------------------------
    template <typename T>
    constexpr const T& min(const T& a, const T& b) noexcept { return (a < b) ? a : b; }

    template <typename T>
    constexpr const T& max(const T& a, const T& b) noexcept { return (a > b) ? a : b; }

    template <typename T>
    constexpr const T& clamp(const T& v, const T& lo, const T& hi) noexcept {
        return (v < lo) ? lo : ((v > hi) ? hi : v);
    }

//----------------------------------------------------
// Sign & Absolute value
//----------------------------------------------------
    template <typename T>
    constexpr int sign(const T& x) noexcept {
        return (T(0) < x) - (x < T(0));
    }

// Generic abs wrapper (works for multiprecision & complex)
    template <typename T>
    inline auto Abs(const T& x) noexcept {
        //using boost::multiprecision::abs; // enable ADL
        using math::abs;
        return abs(x);
    }

//----------------------------------------------------
// Trig / Exp / Log / Sqrt wrappers (ADL-safe)
//----------------------------------------------------
    template <typename T>
    inline auto Sqrt(const T& x) noexcept {
        using std::sqrt;
        return sqrt(x);
    }

    template <typename T>
    inline auto Sin(const T& x) noexcept {
        using std::sin;
        return sin(x);
    }

    template <typename T>
    inline auto Cos(const T& x) noexcept {
        using std::cos;
        return cos(x);
    }

    template <typename T>
    inline auto Tan(const T& x) noexcept {
        using std::tan;
        return tan(x);
    }

    template <typename T>
    inline auto Exp(const T& x) noexcept {
        using std::exp;
        return exp(x);
    }

    template <typename T>
    inline auto Log(const T& x) noexcept {
        using std::log;
        return log(x);
    }

    template <typename Base, typename Exponent>
    inline auto Pow(const Base& a, const Exponent& b) noexcept {
        using std::pow;           // find built-in pow for double/long double
        using boost::multiprecision::pow; // find Boost version for multiprecision
        return pow(a,b);          // unqualified call triggers ADL
    }

    template<typename T>
    inline std::complex<T> PowInt(const std::complex<T>& z, int n) noexcept {
        std::complex<T> result(1,0);
        std::complex<T> base = z;
        int exp = n >= 0 ? n : -n;
        while(exp) {
            if(exp % 2) result *= base;
            base *= base;
            exp /= 2;
        }
        if(n < 0) result = std::complex<T>(1,0)/result;
        return result;
    }

    template<size_t N>
    inline teuk::Real determinant(const std::array<std::array<teuk::Real, N>, N>& M) {
        if constexpr (N == 2) {
            return M[0][0]*M[1][1] - M[0][1]*M[1][0];
        } else if constexpr (N == 3) {
            return M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
                   - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
                   + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
        } else {
            static_assert(N <= 3, "Only 2x2 or 3x3 determinants implemented.");
        }
    }

} // namespace math

// ======================================================
//  Optional: redefine macros to call safe functions
// ======================================================
// #undef SQR
// #define SQR(x) (math::sqr(x))
// #undef CUBE
// #define CUBE(x) (math::cube(x))

#endif // GHZ_NUMERIC_MATHMACROS_HPP


