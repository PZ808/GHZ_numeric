//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_SPINCOEFFSNP_HPP
#define GHZ_NUMERIC_SPINCOEFFSNP_HPP

#pragma once
#include <complex>
#include <map>
#include <string>
#include "TeukTypes.hpp"

using Complex = teuk::Complex;


enum class SpinCoeffType {
    kappa, sigma, lambda, nu,
    rho, mu, tau, pi,
    epsilon, gamma, beta, alpha
};

// NP spin coefficient container
class SpinCoefficients {
private:
    std::map<SpinCoeffType, Complex> coeffs;

public:
    SpinCoefficients() = default;

    void set(SpinCoeffType type, Complex value);
    [[nodiscard]] Complex get(SpinCoeffType type) const;

    // Optional: print all coefficients
    std::string toString() const;
};

#endif //GHZ_NUMERIC_SPINCOEFFSNP_HPP
