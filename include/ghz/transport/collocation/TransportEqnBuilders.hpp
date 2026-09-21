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
#include <algorithm>
#include <limits>
#include <string>

namespace ghz::collocation {

    namespace detail {
        inline void check_rhs_size(const std::vector<teuk::Complex>& rhs,
                                   std::size_t n,
                                   const char* where)
        {
            if (rhs.size() != n) {
                throw std::runtime_error(std::string(where) +
                                         ": right-hand side size mismatch");
            }
        }

        inline teuk::Complex rho_at(teuk::Real r, teuk::Real z, teuk::Real a)
        {
            return -teuk::Real(1) /
                   (teuk::Complex(r, -a * z));
        }
    }

    // Generic source form of the seed equation.  Archive-specific code can
    // condition T_ll and call this without coupling the solver to an archive.
    inline TwoDomainEquationSlice
    build_xmmbar_two_domain_eq_slice(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            teuk::Real z,
            teuk::Real a,
            const std::vector<teuk::Complex>& rhs_left,
            const std::vector<teuk::Complex>& rhs_right)
    {
        using Real = teuk::Real;
        using Complex = teuk::Complex;
        const auto& rL = rops_left.r();
        const auto& rR = rops_right.r();
        detail::check_rhs_size(rhs_left, rL.size(),
                               "build_xmmbar_two_domain_eq_slice(left)");
        detail::check_rhs_size(rhs_right, rR.size(),
                               "build_xmmbar_two_domain_eq_slice(right)");

        TwoDomainEquationSlice eq;
        eq.left.a2.assign(rL.size(), Complex(1.0, 0.0));
        eq.left.a1.resize(rL.size());
        eq.left.a0.resize(rL.size());
        eq.left.rhs = rhs_left;
        eq.right.a2.assign(rR.size(), Complex(1.0, 0.0));
        eq.right.a1.resize(rR.size());
        eq.right.a0.resize(rR.size());
        eq.right.rhs = rhs_right;

        auto fill = [&](const std::vector<Real>& r,
                        ScalarEquationSlice& out) {
            for (std::size_t i = 0; i < r.size(); ++i) {
                const Complex rho = detail::rho_at(r[i], z, a);
                const Complex rhob = std::conj(rho);
                const Real P = (rho + rhob).real();
                const Real Q = ((rho - rhob) * (rho - rhob)).real();
                out.a1[i] = Complex(-P, 0.0);
                out.a0[i] = Complex(-Q, 0.0);
            }
        };
        fill(rL, eq.left);
        fill(rR, eq.right);
        return eq;
    }

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

        return build_xmmbar_two_domain_eq_slice(
                rops_left, rops_right, z, a,
                src_left.Tll_on_r_grid(rops_left.r()),
                src_right.Tll_on_r_grid(rops_right.r()));
    }

    // X_nm equation in direct second-order collocation form.  The source is
    // T_lm + N[X_mmbar] and is deliberately supplied by the caller.
    inline TwoDomainEquationSlice
    build_xnm_two_domain_eq_slice(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            teuk::Real z,
            teuk::Real a,
            const std::vector<teuk::Complex>& rhs_left,
            const std::vector<teuk::Complex>& rhs_right)
    {
        using Complex = teuk::Complex;
        const auto& rL = rops_left.r();
        const auto& rR = rops_right.r();
        detail::check_rhs_size(rhs_left, rL.size(),
                               "build_xnm_two_domain_eq_slice(left)");
        detail::check_rhs_size(rhs_right, rR.size(),
                               "build_xnm_two_domain_eq_slice(right)");

        TwoDomainEquationSlice eq;
        auto fill = [&](const std::vector<teuk::Real>& r,
                        const std::vector<Complex>& rhs,
                        ScalarEquationSlice& out) {
            out.a2.resize(r.size(), Complex(0.5, 0.0));
            out.a1.resize(r.size());
            out.a0.resize(r.size());
            out.rhs = rhs;

            for (std::size_t i = 0; i < r.size(); ++i) {
                const Complex rho = detail::rho_at(r[i], z, a);
                const Complex rhob = std::conj(rho);
                const Complex s = rho + rhob;
                const Complex sp = rho * rho + rhob * rhob;
                const Complex spp = Complex(2.0, 0.0) *
                                    (rho * rho * rho + rhob * rhob * rhob);
                const Complex A = rho * s;
                const Complex Ap = rho * rho * s + rho * sp;
                const Complex App = Complex(2.0, 0.0) * rho * rho * rho * s
                                    + Complex(2.0, 0.0) * rho * rho * sp
                                    + rho * spp;

                // rho/(2s) d_r[s^2 d_r(u/(rho s))] = source.
                out.a1[i] = -rho;
                out.a0[i] = (Ap / A) * (Ap / A)
                            - App / (Complex(2.0, 0.0) * A)
                            - (sp / s) * (Ap / A);
            }
        };
        fill(rL, rhs_left, eq.left);
        fill(rR, rhs_right, eq.right);
        return eq;
    }

    // X_nn is first order after expanding the Held form:
    // 1/2 (rho+rhobar)^2 d_r[X_nn/(rho+rhobar)] = source.
    inline TwoDomainEquationSlice
    build_xnn_two_domain_eq_slice(
            const ghz::numeric::PhysicalChebRadialOps& rops_left,
            const ghz::numeric::PhysicalChebRadialOps& rops_right,
            teuk::Real z,
            teuk::Real a,
            const std::vector<teuk::Complex>& rhs_left,
            const std::vector<teuk::Complex>& rhs_right)
    {
        using Complex = teuk::Complex;
        const auto& rL = rops_left.r();
        const auto& rR = rops_right.r();
        detail::check_rhs_size(rhs_left, rL.size(),
                               "build_xnn_two_domain_eq_slice(left)");
        detail::check_rhs_size(rhs_right, rR.size(),
                               "build_xnn_two_domain_eq_slice(right)");

        TwoDomainEquationSlice eq;
        auto fill = [&](const std::vector<teuk::Real>& r,
                        const std::vector<Complex>& rhs,
                        ScalarEquationSlice& out) {
            out.a2.assign(r.size(), Complex(0.0, 0.0));
            out.a1.resize(r.size());
            out.a0.resize(r.size());
            out.rhs = rhs;
            for (std::size_t i = 0; i < r.size(); ++i) {
                const Complex rho = detail::rho_at(r[i], z, a);
                const Complex rhob = std::conj(rho);
                out.a1[i] = Complex(0.5, 0.0) * (rho + rhob);
                out.a0[i] = -Complex(0.5, 0.0) *
                            (rho * rho + rhob * rhob);
            }
        };
        fill(rL, rhs_left, eq.left);
        fill(rR, rhs_right, eq.right);
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
