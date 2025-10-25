//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_SPINCOEFFSNP_HPP
#define GHZ_NUMERIC_SPINCOEFFSNP_HPP

#pragma once
#include <complex>
#include <map>
#include <string>

using Complex = std::complex<double>;

enum class SpinCoeffType {
    kappa, sigma, lambda, nu,
    rho, mu, tau, pi,
    epsilon, gamma, beta, alpha
};

// GHP Spin Coefficient container
class SpinCoefficients {
private:
    std::map<SpinCoeffType, Complex> coeffs;

public:
    SpinCoefficients() = default;

    void set(SpinCoeffType type, Complex value);
    Complex get(SpinCoeffType type) const;

    // Optional: print all coefficients
    std::string toString() const;
};

#endif //GHZ_NUMERIC_SPINCOEFFSNP_HPP
