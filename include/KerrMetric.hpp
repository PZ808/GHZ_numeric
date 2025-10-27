//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRIC_HPP
#define GHZ_NUMERIC_KERRMETRIC_HPP

#pragma once
#include <array>
#include "Metric.hpp"
#include "VectorsGHZ.hpp"
#include "KerrParams.hpp"


class KerrMetric : public Metric {
private:
    KerrParams params;
    Real k2_, lambda_, mu_, alpha_; // conformal params


public:
    explicit KerrMetric(const KerrParams& p);
    virtual ~KerrMetric() = default;

    Real M() const;
    Real a() const;

    Real r_plus() const;
    Real r_minus() const;

    Real Om_plus() const;
    Real Om_minus() const;

    Real kappa_plus() const;
    Real kappa_minus() const;

    Real k2_C() const;
    Real mu_C() const;
    Real lambda_C() const;
    Real alpha_C () const;

    virtual Real Sigma(double r, double theta) const;
    virtual Real Delta(double r) const;
    virtual Real Lambda(double r, double theta) const;

    Real Sigma_z(Real r, Real z) const;

    Real Lambda_z(Real r, Real z) const;
};

#endif //GHZ_NUMERIC_KERRMETRIC_HPP
