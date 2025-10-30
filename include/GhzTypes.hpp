//
// Created by Peter Zimmerman on 28.10.25.
//

#ifndef GHZ_NUMERIC_GHZTYPES_HPP
#define GHZ_NUMERIC_GHZTYPES_HPP

#pragma once
/**
 * @file GhzTypes.hpp
 * @brief Type and precision definitions for the GHZ_numeric project,
 *        supporting Boost.Multiprecision and runtime precision control.
 */

#include <complex>
#include <limits>
#include <numbers>
#include <type_traits>

// --- Boost Multiprecision ---
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>

namespace teuk {

/**
 * @brief Enumeration for selecting floating-point precision level.
 */
    enum class PrecisionLevel : int {
        Float       = 0,
        Double      = 1,
        LongDouble  = 2,
        Quad        = 3,   ///< 128-bit
        Arbitrary   = 4    ///< arbitrary bit precision (default 256 bits)
    };

/**
 * @brief Compile-time default precision selection.
 * Example: -DTEUK_PRECISION=3 for 128-bit
 */
#ifndef TEUK_PRECISION
#  define TEUK_PRECISION 1
#endif

// ----------------------------------------------------------------------------
// Precision Traits
// ----------------------------------------------------------------------------

    template <PrecisionLevel> struct PrecisionTraits;

// --- float ---
    template <>
    struct PrecisionTraits<PrecisionLevel::Float> {
        using Real     = float;
        using Complex  = std::complex<Real>;
        static constexpr int bits = std::numeric_limits<Real>::digits;
    };

// --- double ---
    template <>
    struct PrecisionTraits<PrecisionLevel::Double> {
        using Real     = double;
        using Complex  = std::complex<Real>;
        static constexpr int bits = std::numeric_limits<Real>::digits;
    };

// --- long double ---
    template <>
    struct PrecisionTraits<PrecisionLevel::LongDouble> {
        using Real     = long double;
        using Complex  = std::complex<Real>;
        static constexpr int bits = std::numeric_limits<Real>::digits;
    };

// --- 128-bit float ---
    template <>
    struct PrecisionTraits<PrecisionLevel::Quad> {
        using Real     = boost::multiprecision::cpp_bin_float_quad;
        using Complex  = std::complex<Real>;
        static constexpr int bits = 113;
    };

// --- Arbitrary (default 256 bits) ---
    template <>
    struct PrecisionTraits<PrecisionLevel::Arbitrary> {

        using Real = boost::multiprecision::number<
                boost::multiprecision::cpp_bin_float<256>,
                boost::multiprecision::et_on>;

        using Complex = std::complex<Real>;
        static constexpr int bits = 256;
    };

// ----------------------------------------------------------------------------
// Runtime precision selection (optional)
// ----------------------------------------------------------------------------

    inline PrecisionLevel active_precision_level =
            static_cast<PrecisionLevel>(TEUK_PRECISION);

    inline void set_precision(PrecisionLevel level) {
        active_precision_level = level;
    }

/**
 * @brief Get the number of bits of precision used by current Real type.
 */
    inline int get_precision_bits() {
        switch (active_precision_level) {
            case PrecisionLevel::Float:       return PrecisionTraits<PrecisionLevel::Float>::bits;
            case PrecisionLevel::Double:      return PrecisionTraits<PrecisionLevel::Double>::bits;
            case PrecisionLevel::LongDouble:  return PrecisionTraits<PrecisionLevel::LongDouble>::bits;
            case PrecisionLevel::Quad:        return PrecisionTraits<PrecisionLevel::Quad>::bits;
            case PrecisionLevel::Arbitrary:   return PrecisionTraits<PrecisionLevel::Arbitrary>::bits;
            default:                          return 53;
        }
    }

// ----------------------------------------------------------------------------
// Type aliases for active precision
// ----------------------------------------------------------------------------

    using ActivePrecision = PrecisionTraits<
            static_cast<PrecisionLevel>(TEUK_PRECISION)>;

    using Real     = typename ActivePrecision::Real;
    using Complex  = typename ActivePrecision::Complex;

    inline constexpr Complex I{0, 1};

// ----------------------------------------------------------------------------
// Utilities
// ----------------------------------------------------------------------------

    inline const Complex CNAN{
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
    };

    template <typename T>
    inline bool Cisnan(const std::complex<T>& z) noexcept {
        using std::isnan;
        return isnan(z.real()) || isnan(z.imag());
    }

    template <typename T>
    inline T Cmaxpart(const std::complex<T>& z) noexcept {
        using std::fabs;
        using std::fmax;
        return fmax(fabs(z.real()), fabs(z.imag()));
    }

// ----------------------------------------------------------------------------
// Generic math wrappers (Boost compatible)
// ----------------------------------------------------------------------------

    template <typename T> inline T Fabs(T x) noexcept { using std::fabs; return fabs(x); }
    template <typename T> inline T Sqrt(T x) noexcept { using std::sqrt; return sqrt(x); }
    template <typename T> inline T Exp(T x)  noexcept { using std::exp;  return exp(x);  }
    template <typename T> inline T Sin(T x)  noexcept { using std::sin;  return sin(x);  }
    template <typename T> inline T Cos(T x)  noexcept { using std::cos;  return cos(x);  }
    template <typename T> inline T Tan(T x)  noexcept { using std::tan;  return tan(x);  }
    template <typename T> inline T Pow(T a, T b) noexcept { using std::pow; return pow(a,b); }

    template <typename T>
    inline std::complex<T> Cexp(const std::complex<T>& z) noexcept {
        using boost::multiprecision::exp;
        return exp(z);
    }

    template <typename T>
    inline std::complex<T> Cpow(const std::complex<T>& a, const std::complex<T>& b) noexcept {
        using boost::multiprecision::pow;
        return pow(a, b);
    }

    template <typename T>
    inline std::complex<T> Csin(const std::complex<T>& z) noexcept {
        using boost::multiprecision::sin;
        return sin(z);
    }

    template <typename T>
    inline std::complex<T> Ccos(const std::complex<T>& z) noexcept {
        using boost::multiprecision::cos;
        return cos(z);
    }

// ----------------------------------------------------------------------------
// Stream operators for multiprecision complex types
// ----------------------------------------------------------------------------

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, const std::complex<T>& z) {
        os << "(" << std::setprecision(20)
           << z.real() << (z.imag() >= 0 ? "+" : "")
           << z.imag() << "i)";
        return os;
    }

} // namespace teuk

namespace teuk::literals {


    //constexpr Real operator""_r(long double x) noexcept { return Real(x); } // real literals
    inline teuk::Real operator ""_r(long double x) noexcept {
        return teuk::Real(x);  // Works for double, long double, Boost multiprecision
    }

    inline teuk::Real operator ""_r(unsigned long long x) noexcept {
        return teuk::Real(x);  // integer literals
    }

// complex literals
// Imaginary literal for teuk::Complex
    inline teuk::Complex operator ""_i(long double x) noexcept {
        return teuk::Complex(0, teuk::Real(x)); // Imaginary part as teuk::Real
    }

    inline teuk::Complex operator ""_i(unsigned long long x) noexcept {
        return teuk::Complex(0, teuk::Real(x)); // integer imaginary literals
    }

    // constexpr Complex operator""_i(long double x) noexcept { return Complex(0, static_cast<Real>(x)); }
// ----------------------
//  full complex literal a + bi
// ----------------------
    inline teuk::Complex operator ""_c(long double x) noexcept {
        // allows writing a + b_c for a purely complex number
        return teuk::Complex(0, teuk::Real(x));
    }
} // namespace teuk::literals

namespace teuk::eigenTypes {
    using Real = double;
    using Complex = std::complex<double>;

    // Dynamic matrix with Real entries
    using MatrixR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    // Dynamics matrix with Complex entries
    using MatrixC = Eigen::Matrix<Complex , Eigen::Dynamic, Eigen::Dynamic>;
    // Dynamic column vector with Complex entries
    using VectorR = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    // Dynamic column vector with Complex entries
    using VectorC = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
}
// namespace teuk::eigenTypes


#endif // GHZ_NUMERIC_GHZTYPES_HPP