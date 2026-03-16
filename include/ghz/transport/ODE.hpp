//
// Created by Peter Zimmerman on 01.12.25.
//
#ifndef GHZ_NUMERIC_ODE_HPP
#define GHZ_NUMERIC_ODE_HPP

#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/core/VectorsGHZ.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <complex>
#include <vector>
#include <functional>
#include <cassert>

using Real = teuk::Real;
using Cplx = teuk::Complex;

namespace ode {
    // State vector type for coupled ODE system at each point z=z_*
    //using StateVec = teuk::CVectorN<4, Real>;   // 3 components: x_mmbar, x_nm, x_nn
    // Here we have to abondon the GHP scalar type and separate into Re and Im parts
    using StateVec = std::vector<Real>; //  components: x_mmbar, Re(x_nm), Im(x_nm), x_nn
    // LHSFunc: computes L(y,r,z_*),
    using LHSFunc    = std::function<StateVec(const StateVec&, Real, Real)>;
    // SourceFunc: computes S(r,y, z_*)
    using SourceFunc = std::function<StateVec(Real, Real)>;

// ==========================================================
// Transport ODE System for dy/dr = LHS(y,r,z) + Source(r,z)
// with full 2D (r,z) integration support.
// ==========================================================
//
///  example usage to solve system dy/dr = L(y,r,z) + S(r,z) of 3 coupled ODEs
//
// std::vector<Real> r_grid = /* r grid */;
// std::vector<Real> z_grid = /* z grid [-1,1] */;

// size_t N = 3; //  e.g. \vec{y} = (x_mmbar, x_nm, x_nn)

// TransportODESystem::StateVec y0(N, 0.0); // initial conditions at r_min
//
//
//auto LHS = [&](const TransportODESystem::StateVec& y, Real r, Real z) -> TransportODESystem::StateVec {
//    TransportODESystem::StateVec dy(N);
//    Cplx rho = -1.0/(r - Cplx(0.0, a*z));
//    // Example: fill dy[i] with LHS expressions
//    dy[0] = (rho + std::conj(rho))*y[1] + (rho - std::conj(rho))*(rho - std::conj(rho))*y[0]; // x_mmbar
//    dy[1] = ...; // x_nm
//    dy[2] = ...; // x_nn
//    return dy;
//};
//
//auto Source = [&](Real r, Real z) -> TransportODESystem::StateVec {
//    TransportODESystem::StateVec s(N);
//    s[0] = T_ll(r,z);
//    s[1] = T_lm(r,z);
//    s[2] = T_ln(r,z);
//    return s;
//};
//
//CoupledODE::TransportODESystem ode{LHS, Source};
//
//auto solution = TransportODESystem::solve_2D(r_grid, z_grid, y0, ode);
// solution[iz][ir][component] gives x_component at (r,z)
////
    struct TransportODESystem {
        LHSFunc LHS;
        SourceFunc Source;

        // functor operator for odeint
        void operator()(const StateVec& y, StateVec& dydr, Real r, Real z) const
        {
            StateVec lhs_eval = LHS(y, r, z);
            StateVec src_eval = Source(r, z);
            assert(dydr.size() == lhs_eval.size() && dydr.size() == src_eval.size());
            for (size_t i=0; i<dydr.size(); ++i)
                dydr[i] = lhs_eval[i] + src_eval[i];
        }
    }; // struct TransportODESystem

    // Adaptive solver for coupled system over r for a single z
    std::vector<StateVec> solve_single_z(
            const std::vector<Real>& r_grid, // nodes in r or (sigma) or the chebyshev grid relative to these values
            const StateVec& y0,              // initial conditions at r_min
            const TransportODESystem& ode,   // ODE system
            Real z,                          // fixed z value
            Real abs_tol=1e-9,             // solver abs tolerance
            Real rel_tol=1e-9,             // solver rel tolerance
            Real initial_step=1e-3);      // initial step size

} // namespace ode


#endif //GHZ_NUMERIC_ODE_HPP
