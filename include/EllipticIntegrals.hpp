//
// Created by Peter Zimmerman on 02.11.25.
//

#ifndef GHZ_NUMERIC_ELLIPTICINTEGRALS_HPP
#define GHZ_NUMERIC_ELLIPTICINTEGRALS_HPP

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "GhzTypes.hpp"

using Real = teuk::Real;


class EllipticIntegrals {
public:
    // Constructor with no fixed parameters
    EllipticIntegrals() {}

    // Function to compute the complete elliptic integral of the first kind F(k)
    static Real computeFirstKind(const Real& k) {
        return boost::math::ellint_1(k);
    }

    // Function to compute the complete elliptic integral of the second kind E(k)
    static Real computeSecondKind(const Real& k) {
        return boost::math::ellint_2(k);
    }

    // Function to compute the complete elliptic integral of the third kind Pi(n, k)
    static Real computeThirdKind(const Real& n, const Real& k) {
        return boost::math::ellint_3(n, k);
    }
};

#endif //GHZ_NUMERIC_ELLIPTICINTEGRALS_HPP
