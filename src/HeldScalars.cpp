//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/HeldScalars.hpp"

#pragma once
#include "SpinCoefficients.hpp"
#include "Tetrad.hpp"

#include "../include/HeldScalars.hpp"
#include <iostream>

void HeldScalars::fromSpinCoeffs(const SpinCoefficients& sc) {
    using SCT = SpinCoeffType;

    // Unprimed
    rho = sc.get(SCT::rho);
    tau = sc.get(SCT::tau);
    epsilon = sc.get(SCT::epsilon);
    kappa = sc.get(SCT::kappa);
    sigma = sc.get(SCT::sigma);

    // Primed (Held conjugate)
    rhoPrime = -std::conj(sc.get(SCT::mu));
    tauPrime = std::conj(sc.get(SCT::pi));
    epsilonPrime = std::conj(sc.get(SCT::gamma));
    kappaPrime = -std::conj(sc.get(SCT::nu));
    sigmaPrime = -std::conj(sc.get(SCT::lambda));
}

void HeldScalars::print() const {
    auto printC = [](std::string name, Complex c) {
        std::cout << name << " = " << c << "\n";
    };

    printC("rho", rho);
    printC("tau", tau);
    printC("epsilon", epsilon);
    printC("kappa", kappa);
    printC("sigma", sigma);

    printC("rho'", rhoPrime);
    printC("tau'", tauPrime);
    printC("epsilon'", epsilonPrime);
    printC("kappa'", kappaPrime);
    printC("sigma'", sigmaPrime);
}
