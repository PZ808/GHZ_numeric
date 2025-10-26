//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/SpinCoeffsNP.hpp"

void SpinCoefficients::set(SpinCoeffType type, Complex value) {
    coeffs[type] = value;
}

Complex SpinCoefficients::get(SpinCoeffType type) const {
    auto it = coeffs.find(type);
    if(it != coeffs.end()) return it->second;
    return 0.0;
}