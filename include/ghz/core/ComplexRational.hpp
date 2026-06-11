//
// Created by Peter Zimmerman on 24.04.26.
//

#ifndef GHZ_NUMERIC_COMPLEXRATIONAL_HPP
#define GHZ_NUMERIC_COMPLEXRATIONAL_HPP

#pragma once

#include "Floating.hpp"
#include "GhzTypes.hpp"

#include <complex>
#include <span>
#include <vector>

namespace teuk::cx {

    template <class Real>
    inline Real abs2(const std::complex<Real>& z) {
        return z.real() * z.real() + z.imag() * z.imag();
    }

    template <class Real>
    inline Real re_conj(const std::complex<Real>& a, const std::complex<Real>& b) {
        return (a * std::conj(b)).real();
    }

    template <class Real>
    inline Real im_conj(const std::complex<Real>& a, const std::complex<Real>& b) {
        return (a * std::conj(b)).imag();
    }

    template <class Real>
    inline std::complex<Real> inverse_checked(
            const std::complex<Real>& z,
            const char* where,
            const Real& multiplier = teuk::floating::default_tolerance_multiplier<Real>(),
            const Real& floor = teuk::floating::default_scale_floor<Real>())
    {
        const Real denom = teuk::floating::denominator_checked(
                abs2(z), where, Real(1), multiplier, floor);
        return std::conj(z) / denom;
    }

    template <class Real>
    inline std::complex<Real> divide_checked(
            const std::complex<Real>& num,
            const std::complex<Real>& den,
            const char* where,
            const Real& multiplier = teuk::floating::default_tolerance_multiplier<Real>(),
            const Real& floor = teuk::floating::default_scale_floor<Real>())
    {
        return num * inverse_checked(den, where, multiplier, floor);
    }

    template <class Real>
    inline std::complex<Real> divide_by_real_checked(
            const std::complex<Real>& num,
            const Real& den,
            const char* where,
            const Real& multiplier = teuk::floating::default_tolerance_multiplier<Real>(),
            const Real& floor = teuk::floating::default_scale_floor<Real>())
    {
        const Real d = teuk::floating::denominator_checked(
                den, where, Real(1), multiplier, floor);
        return num / d;
    }

    template <class Real>
    inline Real real_ratio_checked(
            const std::complex<Real>& num,
            const std::complex<Real>& den,
            const char* where,
            const Real& multiplier = teuk::floating::default_tolerance_multiplier<Real>(),
            const Real& floor = teuk::floating::default_scale_floor<Real>())
    {
        const auto q = divide_checked(num, den, where, multiplier, floor);
        return teuk::floating::real_part_checked(q, where, multiplier, floor);
    }

    template <class Real>
    inline std::complex<Real> project_to_real(
            const std::complex<Real>& z,
            const char* where,
            const Real& multiplier = teuk::floating::default_tolerance_multiplier<Real>(),
            const Real& floor = teuk::floating::default_scale_floor<Real>())
    {
        return std::complex<Real>(
                teuk::floating::real_part_checked(z, where, multiplier, floor),
                Real(0));
    }

    template <class Real>
    inline std::complex<Real> project_to_imag(
            const std::complex<Real>& z,
            const char* where,
            const Real& multiplier = teuk::floating::default_tolerance_multiplier<Real>(),
            const Real& floor = teuk::floating::default_scale_floor<Real>())
    {
        using std::abs;
        const Real scale = std::max(abs(z.real()), std::max(abs(z.imag()), floor));
        const Real tol = teuk::floating::scaled_tolerance(scale, multiplier, floor);
        if (abs(z.real()) > tol) {
            throw std::runtime_error(std::string(where) + ": expected nearly imaginary quantity.");
        }
        return std::complex<Real>(Real(0), z.imag());
    }

    template <class Real>
    inline bool nearly_conjugate(
            const std::complex<Real>& a,
            const std::complex<Real>& b,
            const Real& multiplier = teuk::floating::default_tolerance_multiplier<Real>(),
            const Real& floor = teuk::floating::default_scale_floor<Real>())
    {
        return teuk::floating::nearly_equal(a, std::conj(b), multiplier, floor);
    }

    template <class Real>
    inline std::complex<Real> mobius_checked(
            const std::complex<Real>& a,
            const std::complex<Real>& b,
            const std::complex<Real>& c,
            const std::complex<Real>& d,
            const std::complex<Real>& z,
            const char* where,
            const Real& multiplier = teuk::floating::default_tolerance_multiplier<Real>(),
            const Real& floor = teuk::floating::default_scale_floor<Real>())
    {
        return divide_checked(a * z + b, c * z + d, where, multiplier, floor);
    }

    template <class Complex>
    Complex horner_eval(std::span<const Complex> coeffs, const Complex& z) {
        Complex out = Complex(0);
        for (const auto& c : coeffs) out = out * z + c;
        return out;
    }

    template <class Complex>
    Complex rational_horner_checked(
            std::span<const Complex> num_coeffs,
            std::span<const Complex> den_coeffs,
            const Complex& z,
            const char* where)
    {
        const Complex num = horner_eval(num_coeffs, z);
        const Complex den = horner_eval(den_coeffs, z);
        return divide_checked(num, den, where);
    }

    template <class Real>
    inline std::complex<Real> rho_plus_rhob(const std::complex<Real>& rho) {
        return rho + std::conj(rho);
    }

    template <class Real>
    inline std::complex<Real> rho_minus_rhob(const std::complex<Real>& rho) {
        return rho - std::conj(rho);
    }

    template <class Real>
    inline Real real_part_of_rho_plus_rhob_checked(
            const std::complex<Real>& rho,
            const char* where)
    {
        return teuk::floating::denominator_checked(rho + std::conj(rho), where);
    }


    template <class Real>
    inline Real imag_rho_minus_rhob_checked(
            const std::complex<Real>& rho,
            const char* where)
    {
        const std::complex<Real> diff = rho - std::conj(rho);
        return teuk::floating::imag_part_checked(diff, where);
    }

    template <class Real>
    inline Real checked_imag_rho_minus_rhob(
            const std::complex<Real>& rho,
            const char* where)
    {
        const std::complex<Real> diff = rho - std::conj(rho);
        using std::abs;
        const Real tol = teuk::floating::scaled_tolerance(
                std::max(abs(diff.real()), std::max(abs(diff.imag()), Real(1)))
        );
        if (abs(diff.real()) > tol) {
            throw std::runtime_error(std::string(where) + ": expected nearly imaginary rho-rhob.");
        }
        return diff.imag();
    }

    template <class Real>
    inline std::complex<Real> reciprocal_rho_sum_checked(
            const std::complex<Real>& rho,
            const char* where)
    {
        const Real denom = teuk::floating::denominator_checked(
                rho + std::conj(rho), where);
        return std::complex<Real>(Real(1) / denom, Real(0));
    }

} // namespace teuk::cx

#endif //GHZ_NUMERIC_COMPLEXRATIONAL_HPP
