//
// Created by Peter Zimmerman on 24.04.26.
//
#ifndef GHZ_NUMERIC_FLOATING_HPP
#define GHZ_NUMERIC_FLOATING_HPP

#pragma once
/**
 * @file Floating.hpp
 * @brief Precision-aware floating-point utilities for GHZ_numeric.
 *
 * This header centralizes tolerance handling and floating-point sanity checks.
 * It is designed to work with teuk::Real / teuk::Complex from GhzTypes.hpp,
 * including Boost.Multiprecision real types.
 */

#include "GhzTypes.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#ifndef TEUK_FLOATING_DEFAULT_TOL_MULTIPLIER
#define TEUK_FLOATING_DEFAULT_TOL_MULTIPLIER 100
#endif

#ifndef TEUK_FLOATING_DEFAULT_SCALE_FLOOR
#define TEUK_FLOATING_DEFAULT_SCALE_FLOOR 1
#endif

namespace teuk::floating {

// -----------------------------------------------------------------------------
// Type traits
// -----------------------------------------------------------------------------

    template <class T>
    struct RealType {
        using type = T;
    };

    template <class T>
    struct RealType<std::complex<T>> {
        using type = T;
    };

    template <class T>
    using RealType_t = typename RealType<T>::type;

// -----------------------------------------------------------------------------
// Defaults tied to GHZ types
// -----------------------------------------------------------------------------

    using DefaultReal    = teuk::Real;
    using DefaultComplex = teuk::Complex;

// -----------------------------------------------------------------------------
// Basic helpers
// -----------------------------------------------------------------------------

    template <class Real>
    inline Real machine_epsilon() {
        return std::numeric_limits<Real>::epsilon();
    }

    template <class Real>
    inline Real default_tolerance_multiplier() {
        return static_cast<Real>(TEUK_FLOATING_DEFAULT_TOL_MULTIPLIER);
    }

    template <class Real>
    inline Real default_scale_floor() {
        return static_cast<Real>(TEUK_FLOATING_DEFAULT_SCALE_FLOOR);
    }

    template <class Real>
    inline Real abs_value(const Real& x) {
        using std::abs;
        return abs(x);
    }

    template <class Real>
    inline Real abs_value(const std::complex<Real>& z) {
        using std::abs;
        return abs(z);
    }

    template <class Real>
    inline Real max3(const Real& a, const Real& b, const Real& c) {
        return std::max(a, std::max(b, c));
    }

    template <class Real>
    inline Real scaled_tolerance(
            const Real& scale,
            const Real& multiplier = default_tolerance_multiplier<Real>(),
            const Real& floor = default_scale_floor<Real>())
    {
        return multiplier * machine_epsilon<Real>() * std::max(scale, floor);
    }

// -----------------------------------------------------------------------------
// Generic tolerance builders
// -----------------------------------------------------------------------------

    template <class T>
    inline RealType_t<T> tolerance_for(
            const T& x,
            const RealType_t<T>& multiplier = default_tolerance_multiplier<RealType_t<T>>(),
            const RealType_t<T>& floor = default_scale_floor<RealType_t<T>>())
    {
        using Real = RealType_t<T>;
        return scaled_tolerance<Real>(abs_value(x), multiplier, floor);
    }

    template <class T>
    inline RealType_t<T> pair_tolerance(
            const T& a,
            const T& b,
            const RealType_t<T>& multiplier = default_tolerance_multiplier<RealType_t<T>>(),
            const RealType_t<T>& floor = default_scale_floor<RealType_t<T>>())
    {
        using Real = RealType_t<T>;
        const Real scale = max3<Real>(abs_value(a), abs_value(b), floor);
        return scaled_tolerance<Real>(scale, multiplier, floor);
    }

// -----------------------------------------------------------------------------
// Predicates
// -----------------------------------------------------------------------------

    template <class T>
    inline bool nearly_equal(
            const T& a,
            const T& b,
            const RealType_t<T>& multiplier = default_tolerance_multiplier<RealType_t<T>>(),
            const RealType_t<T>& floor = default_scale_floor<RealType_t<T>>())
    {
        return abs_value(a - b) <= pair_tolerance(a, b, multiplier, floor);
    }

    template <class T>
    inline bool nearly_zero(
            const T& x,
            const RealType_t<T>& reference_scale = default_scale_floor<RealType_t<T>>(),
            const RealType_t<T>& multiplier = default_tolerance_multiplier<RealType_t<T>>(),
            const RealType_t<T>& floor = default_scale_floor<RealType_t<T>>())
    {
        using Real = RealType_t<T>;
        const Real scale = max3<Real>(abs_value(x), reference_scale, floor);
        return abs_value(x) <= scaled_tolerance<Real>(scale, multiplier, floor);
    }

    template <class Real>
    inline bool nearly_real(
            const std::complex<Real>& z,
            const Real& multiplier = default_tolerance_multiplier<Real>(),
            const Real& floor = default_scale_floor<Real>())
    {
        using std::abs;
        const Real scale = max3<Real>(abs(z.real()), abs(z.imag()), floor);
        return abs(z.imag()) <= scaled_tolerance<Real>(scale, multiplier, floor);
    }

// -----------------------------------------------------------------------------
// Checked extraction / denominators
// -----------------------------------------------------------------------------

    template <class Real>
    inline Real real_part_checked(
            const std::complex<Real>& z,
            const char* where,
            const Real& multiplier = default_tolerance_multiplier<Real>(),
            const Real& floor = default_scale_floor<Real>())
    {
        if (!nearly_real(z, multiplier, floor)) {
            std::ostringstream oss;
            oss << where << ": expected a nearly real quantity, got "
                << z.real() << " + i*" << z.imag();
            throw std::runtime_error(oss.str());
        }
        return z.real();
    }

    template <class Real>
    inline Real denominator_checked(
            const Real& x,
            const char* where,
            const Real& reference_scale = default_scale_floor<Real>(),
            const Real& multiplier = default_tolerance_multiplier<Real>(),
            const Real& floor = default_scale_floor<Real>())
    {
        using std::abs;
        const Real scale = max3<Real>(abs(x), reference_scale, floor);
        const Real tol   = scaled_tolerance<Real>(scale, multiplier, floor);

        if (abs(x) <= tol) {
            std::ostringstream oss;
            oss << where << ": denominator too small. value=" << x
                << ", tol=" << tol;
            throw std::runtime_error(oss.str());
        }

        return x;
    }

    template <class Real>
    inline Real denominator_checked(
            const std::complex<Real>& z,
            const char* where,
            const Real& reference_scale = default_scale_floor<Real>(),
            const Real& multiplier = default_tolerance_multiplier<Real>(),
            const Real& floor = default_scale_floor<Real>())
    {
        const Real x = real_part_checked(z, where, multiplier, floor);
        return denominator_checked<Real>(x, where, reference_scale, multiplier, floor);
    }

    template <class Real>
    inline Real imag_part_checked(
            const std::complex<Real>& z,
            const char* where,
            const Real& multiplier = default_tolerance_multiplier<Real>(),
            const Real& floor = default_scale_floor<Real>())
    {
        using std::abs;
        const Real scale = std::max(abs(z.real()), std::max(abs(z.imag()), floor));
        const Real tol   = scaled_tolerance<Real>(scale, multiplier, floor);

        if (abs(z.real()) > tol) {
            std::ostringstream oss;
            oss << where << ": expected a nearly imaginary quantity, got "
                << z.real() << " + i*" << z.imag();
            throw std::runtime_error(oss.str());
        }

        return z.imag();
    }
// -----------------------------------------------------------------------------
// Throwing assertion helpers
// -----------------------------------------------------------------------------

    template <class T>
    inline void require_nearly_equal(
            const T& a,
            const T& b,
            const char* where,
            const RealType_t<T>& multiplier = default_tolerance_multiplier<RealType_t<T>>(),
            const RealType_t<T>& floor = default_scale_floor<RealType_t<T>>())
    {
        if (!nearly_equal(a, b, multiplier, floor)) {
            const auto diff = abs_value(a - b);
            const auto tol  = pair_tolerance(a, b, multiplier, floor);

            std::ostringstream oss;
            oss << where << ": values are not nearly equal. "
                << "|a-b|=" << diff << ", tol=" << tol;
            throw std::runtime_error(oss.str());
        }
    }

    template <class T>
    inline void require_nearly_zero(
            const T& x,
            const char* where,
            const RealType_t<T>& reference_scale = default_scale_floor<RealType_t<T>>(),
            const RealType_t<T>& multiplier = default_tolerance_multiplier<RealType_t<T>>(),
            const RealType_t<T>& floor = default_scale_floor<RealType_t<T>>())
    {
        using Real = RealType_t<T>;
        const Real scale = max3<Real>(abs_value(x), reference_scale, floor);
        const Real tol   = scaled_tolerance<Real>(scale, multiplier, floor);

        if (abs_value(x) > tol) {
            std::ostringstream oss;
            oss << where << ": quantity is not nearly zero. "
                << "|x|=" << abs_value(x) << ", tol=" << tol;
            throw std::runtime_error(oss.str());
        }
    }

// -----------------------------------------------------------------------------
// Optional policy object
// -----------------------------------------------------------------------------

    template <class Real = DefaultReal>
    struct TolerancePolicy {
        Real multiplier = default_tolerance_multiplier<Real>();
        Real floor      = default_scale_floor<Real>();
    };

} // namespace teuk::floating

#endif //GHZ_NUMERIC_FLOATING_HPP