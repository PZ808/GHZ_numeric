//
// Created by Peter Zimmerman on 28.10.25.
//

#ifndef GHZ_NUMERIC_EIGENTEUKNUMTRAITS_HPP
#define GHZ_NUMERIC_EIGENTEUKNUMTRAITS_HPP

// EigenTeukNumTraits.hpp
#pragma once
//#include <Eigen/Core>
#include "GhzTypes.hpp" // defines teuk::Real and teuk::Complex
#include <complex>

namespace Eigen {

// For multiprecision real
    template<>
    struct NumTraits<teuk::Real> : GenericNumTraits<teuk::Real>
    {
        using RealType = teuk::Real;
        using Literal  = teuk::Real;
        using NonInteger = teuk::Real;
        using Nested  = teuk::Real;

        enum {
            IsComplex = 0,
            IsInteger = 0,
            IsSigned  = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost  = 2,
            MulCost  = 2
        };
    };

// For multiprecision complex
    template<>
    struct NumTraits<teuk::Complex> : GenericNumTraits<teuk::Complex>
    {
        using RealType = teuk::Real;
        using Literal  = teuk::Complex;
        using NonInteger = teuk::Complex;
        using Nested  = teuk::Complex;

        enum {
            IsComplex = 1,
            IsInteger = 0,
            IsSigned  = 1,
            RequireInitialization = 1,
            ReadCost = 2,
            AddCost  = 4,
            MulCost  = 4
        };
    };

} // namespace Eigen


namespace Eigen {

#if TEUK_PRECISION > 2

    #include <boost/multiprecision/cpp_bin_float.hpp>
    #include <complex>

    namespace Eigen {

    template<>
    struct NumTraits<teuk::Real> : GenericNumTraits<teuk::Real>
    {
        using Real = teuk::Real;
        using Literal = Real;
        using NonInteger = Real;
        using Nested = Real;

        enum {
            IsComplex = 0,
            IsInteger = 0,
            IsSigned  = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost  = 2,
            MulCost  = 2
        };
    };

    template<>
    struct NumTraits<teuk::Complex> : GenericNumTraits<teuk::Complex>
    {
        using Real = teuk::Real;
        using Literal = teuk::Complex;
        using NonInteger = teuk::Complex;
        using Nested = teuk::Complex;

        enum {
            IsComplex = 1,
            IsInteger = 0,
            IsSigned  = 1,
            RequireInitialization = 1,
            ReadCost = 2,
            AddCost  = 4,
            MulCost  = 4
        };
    };


#endif // TEUK_PRECISION > 2
} // namespace Eigen

#endif //GHZ_NUMERIC_EIGENTEUKNUMTRAITS_HPP
