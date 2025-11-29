//
// Created by Peter Zimmerman on 25.10.25.
//

#include <utility>

#include "../include/SpinCoeffsNP.hpp"

void SpinCoefficients::set(SpinCoeffType type, Complex value) {
    coeffs[type] = std::move(value);
}

Complex SpinCoefficients::get(SpinCoeffType type) const {
    auto it = coeffs.find(type);
    if (it != coeffs.end()) return it->second;
    return teuk::zeroC;
}