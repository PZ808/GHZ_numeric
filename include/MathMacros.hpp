//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_MATHMACROS_HPP
#define GHZ_NUMERIC_MATHMACROS_HPP

// simple_math_macros.hpp
// Small collection of simple math macros and their safer inline constexpr alternatives.
//
// Usage:
//   #include "math_macros.hpp"
//
//   auto a = SQR(x);        // macro (be careful with side-effects)
//   auto b = math::sqr(x);  // safer, constexpr template function
//
// Notes:
// - Macros evaluate their arguments multiple times; prefer the inline constexpr functions below.
// - The macro names are provided for compatibility with C-style code that expects macros.

// Include <cmath> only for std::pow or std::abs when needed.
#include <cmath>
#include <type_traits>

//
// MACROS (simple, but may evaluate arguments multiple times)
// Always parenthesize arguments and the whole expression.
//
#define SQR(x)    ((x) * (x))
#define CUBE(x)   ((x) * (x) * (x))
#define POW2(x)   SQR(x)   // alias
#define POW3(x)   CUBE(x)  // alias

// MIN / MAX macros (beware of double evaluation)
#ifndef MIN
#define MIN(a, b) (( (a) < (b) ) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) (( (a) > (b) ) ? (a) : (b))
#endif

// CLAMP macro: clamps v to [lo, hi]
#ifndef CLAMP
#define CLAMP(v, lo, hi) ( ((v) < (lo)) ? (lo) : (((v) > (hi)) ? (hi) : (v)) )
#endif

// SIGN macro: -1, 0, +1
#ifndef SIGN
#define SIGN(x) ( ((x) > 0) - ((x) < 0) )
#endif

//
// Safer, type-generic, constexpr inline alternatives (preferred in C++)
//
namespace math {

// sqr: x*x
    template<typename T>
    constexpr T sqr(T x) noexcept
    {
        return x * x;
    }

// cube: x*x*x
    template<typename T>
    constexpr T cube(T x) noexcept
    {
        return x * x * x;
    }

// pow_n: integer small exponents (2..4). For general exponents use std::pow (floating).
    template<typename T>
    constexpr T pow2(T x) noexcept { return x * x; }
    template<typename T>
    constexpr T pow3(T x) noexcept { return x * x * x; }
    template<typename T>
    constexpr T pow4(T x) noexcept { return (x * x) * (x * x); }

// min / max
    template<typename T>
    constexpr const T& min(const T& a, const T& b) noexcept { return (a < b) ? a : b; }
    template<typename T>
    constexpr const T& max(const T& a, const T& b) noexcept { return (a > b) ? a : b; }

// clamp
    template<typename T>
    constexpr const T& clamp(const T& v, const T& lo, const T& hi) noexcept
    {
        return (v < lo) ? lo : ((v > hi) ? hi : v);
    }

// sign: -1, 0, +1 (integral and floating)
    template<typename T>
    constexpr int sign(T x) noexcept
    {
        return (T(0) < x) - (x < T(0));
    }

// safe abs wrapper that works for integers and floating point
    template<typename T>
    constexpr T absval(T x) noexcept
    {
        return x < T(0) ? -x : x;
    }

} // namespace math

//
// Optional: if you still want macro names that forward to safer functions,
// uncomment the following lines:
//
// #undef SQR
// #define SQR(x) (math::sqr(x))
//
// This keeps the original macro name but uses the inline function, avoiding
// multiple evaluations. Only do this if macros are required by existing code.
//

#endif //GHZ_NUMERIC_MATHMACROS_HPP
