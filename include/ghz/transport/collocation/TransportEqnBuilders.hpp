//
// Created by Peter Zimmerman on 19.03.26.
//

#ifndef HS1DATA_DATA_TRANSPORTEQNBUILDERS_HPP
#define HS1DATA_DATA_TRANSPORTEQNBUILDERS_HPP

#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/source/ConditionSource.hpp"
#include "ghz/transport/collocation/SpectralSliceSolve.hpp"
#include "ghz/spectral/PhysicalChebRadialOps.hpp"

#include <stdexcept>
#include <vector>
#include <complex>

namespace ghz::collocation {

    /**
     * @brief Build the two-domain collocation equation for X_{m\bar m}
     * on a fixed z-slice, assuming the conditioned sources already have
     * that z-slice active.
     *
     * This corresponds to the old RK builder
     *
     *   u'' = P(r,z) u' + Q(r,z) u + T_ll(r,z),
     *
     * with
     *
     *   P = Re(rho + rhob),
     *   Q = Re((rho - rhob)^2),
     *   rho = -1 / (r - i a z).
     *
     * So in collocation form
     *
     *   a2 u'' + a1 u' + a0 u = rhs
     *
     * we use
     *
     *   a2 = 1,
     *   a1 = -p,
     *   a0 = -q,
     *   rhs = T_ll.
     *
     * @param rops_left   physical Chebyshev radial ops on [r_min, r_p]
     * @param rops_right  physical Chebyshev radial ops on [r_p, r_max]
     * @param z           fixed z-slice value
     * @param a           Kerr spin parameter
     * @param src_left    conditioned source on left patch, already set to z
     * @param src_right   conditioned source on right patch, already set to z
     */
    inline TwoDomainEquationSlice
    build_xmmbar_two_domain_eq_slice_prepared(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            teuk::Real z,
            teuk::Real a,
            const ghz::source::ConditionSource& src_left,
            const ghz::source::ConditionSource& src_right)
    {
        using Real = teuk::Real;
        using Complex = teuk::Complex;

        if (!src_left.has_active_z_slice()) {
            throw std::runtime_error(
                    "build_xmmbar_two_domain_eq_slice_prepared: left source has no active z-slice.");
        }
        if (!src_right.has_active_z_slice()) {
            throw std::runtime_error(
                    "build_xmmbar_two_domain_eq_slice_prepared: right source has no active z-slice.");
        }

        // Strong sanity check: make sure the conditioned source slices match z.
        {
            const Real zl = src_left.z_slice();
            const Real zr = src_right.z_slice();
            const Real tol = Real(100) * std::numeric_limits<Real>::epsilon();

            if (std::abs(zl - z) > tol * std::max(Real(1), std::abs(z))) {
                throw std::runtime_error(
                        "build_xmmbar_two_domain_eq_slice_prepared: left source z-slice does not match requested z.");
            }
            if (std::abs(zr - z) > tol * std::max(Real(1), std::abs(z))) {
                throw std::runtime_error(
                        "build_xmmbar_two_domain_eq_slice_prepared: right source z-slice does not match requested z.");
            }
        }

        const auto& rL = rops_left.r();
        const auto& rR = rops_right.r();

        const std::size_t NL = rL.size();
        const std::size_t NR = rR.size();

        TwoDomainEquationSlice eq;
        eq.left.a2.resize(NL);
        eq.left.a1.resize(NL);
        eq.left.a0.resize(NL);
        eq.left.rhs = src_left.Tll_on_r_grid(rL);

        eq.right.a2.resize(NR);
        eq.right.a1.resize(NR);
        eq.right.a0.resize(NR);
        eq.right.rhs = src_right.Tll_on_r_grid(rR);

        // Left coefficients
        for (std::size_t i = 0; i < NL; ++i) {
            const Real r = rL[i];
            const Complex rho  = -Real(1) / (r - Complex(0.0, a*z));
            const Complex rhob = std::conj(rho);

            const Real P = (rho + rhob).real();
            const Real Q = ((rho - rhob) * (rho - rhob)).real();

            eq.left.a2[i] = Complex(1.0, 0.0);
            eq.left.a1[i] = Complex(-P, 0.0);
            eq.left.a0[i] = Complex(-Q, 0.0);
        }

        // Right coefficients
        for (std::size_t i = 0; i < NR; ++i) {
            const Real r = rR[i];
            const Complex rho  = -Real(1) / (r - Complex(0.0, a * z));
            const Complex rhob = std::conj(rho);

            const Real P = (rho + rhob).real();
            const Real Q = ((rho - rhob) * (rho - rhob)).real();

            eq.right.a2[i] = Complex(1.0, 0.0);
            eq.right.a1[i] = Complex(-P, 0.0);
            eq.right.a0[i] = Complex(-Q, 0.0);
        }

        return eq;
    }

    /**
     * @brief Convenience wrapper: activates the z-slice on both conditioned
     * sources and then builds the two-domain X_{m\bar m} equation.
     *
     * Use the *_prepared version inside a solve loop if you already call
     * set_z_slice(z) outside for efficiency.
     */
    inline TwoDomainEquationSlice
    build_xmmbar_two_domain_eq_slice(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            teuk::Real z,
            teuk::Real a,
            ghz::source::ConditionSource& src_left,
            ghz::source::ConditionSource& src_right)
    {
        src_left.set_z_slice(z);
        src_right.set_z_slice(z);

        return build_xmmbar_two_domain_eq_slice_prepared(
                rops_left, rops_right, z, a, src_left, src_right);
    }

} // namespace ghz::collocation

#endif //HS1DATA_DATA_TRANSPORTEQNBUILDERS_HPP
