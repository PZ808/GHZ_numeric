//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_GHPSCALARS_HPP
#define GHZ_NUMERIC_GHPSCALARS_HPP


struct GHPType {
    int p;
    int q;
    GHPType(int p_in, int q_in) : p(p_in), q(q_in) {}
};


class GHPScalarsTypeD {
    // GHP scalars for type-D spacetime in a Weyl aligned NP frame
    // (shear-free geodesic null congruence)
public:
    // GHP scalars (unprimed and primed)
    Complex rho, rhop, rho_bar, rhop_bar;
    Complex tau, taup, tau_bar, taup_bar;
    Complex Psi2, Psi2_bar;

    GHPType rho_GHPType{1, 1};
    GHPType rho_bar_GHPType{1, 1};

    GHPType rhop_GHPType{-1,-1};
    GHPType rhop_bar_GHPType{-1,-1};

    GHPType tau_GHPType{1, -1};
    GHPType tau_bar_GHPType{-1, 1};

    GHPType taup_GHPType{-1, 1};
    GHPType taup_bar_GHPType{1, -1};

    GHPType Psi2_GHPType{0, 0};
    GHPType Psi2_bar_GHPType{0, 0};

    GHPScalarsTypeD() = default;

    void fromSpinCoeffs(const SpinCoefficients& sc);
};

#endif //GHZ_NUMERIC_GHPSCALARS_HPP
