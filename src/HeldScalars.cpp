//
// Created by Peter Zimmerman on 25.10.25.
//

#pragma once
#include <iostream>
#include "../include/SpinCoeffsNP.hpp"
#include "../include/GHPScalars.hpp"
#include "../include/WeylScalars.hpp"
#include "../include/HeldScalars.hpp"
#include <cmath>


void HeldScalars::set_HeldScalars_from_NP_SpinCoeffs(const SpinCoefficients& sc, const WeylScalars& ws) {
    using SCT = SpinCoeffType;
    using WST = WeylScalarType;
    using std::conj;
    // Unprimed
    Complex rho = sc.get(SCT::rho);
    Complex rho_bar = conj(rho);

    tauH = sc.get(SCT::tau)/(rho*rho_bar);
    tauH_bar = conj(tauH);

    rhopH = -0.5;
    rhopH_bar = -0.5;

    PsiH = ws.get(WST::Psi2)/(rho*rho*rho);

}

void HeldScalars::print() const {
    auto printC = [](std::string name, Complex c) {
        std::cout << name << " = " << c << "\n";
    };

    printC("tauH", tauH);
    printC("rhopH", rhopH);
    printC("PsiH", PsiH);
    printC("OmH", OmH);

}
