//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_WEYLSCALARS_HPP
#define GHZ_NUMERIC_WEYLSCALARS_HPP

#pragma once
#include <complex>
#include <map>
#include "TeukTypes.hpp"

enum class WeylScalarType {
    Psi0, Psi1, Psi2, Psi3, Psi4
};


// GHP Spin Coefficient container
class WeylScalars {
    using Complex = teuk::Complex;

private:
    std::map<WeylScalarType, Complex> weyl_scalars;

public:
    WeylScalars() = default;

    void set(WeylScalarType type, Complex value);
    Complex get(WeylScalarType type) const;

    // Optional: print all coefficients
    std::string toString() const;
};

#endif //GHZ_NUMERIC_WEYLSCALARS_HPP
