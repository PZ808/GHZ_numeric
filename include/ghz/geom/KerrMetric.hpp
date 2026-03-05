//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRIC_HPP
#define GHZ_NUMERIC_KERRMETRIC_HPP

#pragma once
#include <array>
#include "ghz/geom/Metric.hpp"
#include "ghz/core/VectorsGHZ.hpp"
#include "KerrParams.hpp"


class KerrMetric : public Metric {
private:
    KerrParams params;
    Real k2_, lambda_, mu_, alpha_; // conformal params


public:
    explicit KerrMetric(const KerrParams& p);
    virtual ~KerrMetric() = default;

    //
    // background parameters
    //
    Real M() const; // mass
    Real a() const; // spin
    Real r_plus() const; // outer horizon
    Real r_minus() const; // inner horizon
    Real Om_plus() const; // outer horizon freq
    Real Om_minus() const; // inner horizon rotn freq
    Real kappa_plus() const; // outer horizon surf. grav.
    Real kappa_minus() const; // inner horizon surf. grav.
    // conformal paramaters see https://arxiv.org/pdf/1910.13452
    Real k2_C() const;
    Real mu_C() const;
    Real lambda_C() const;
    Real alpha_C () const;
    //
    // metric functions
    //
    virtual Real Sigma(Real r, Real theta) const;
    virtual Real Delta(Real r) const;
    virtual Real Lambda(Real r, Real theta) const;
    Real Sigma_z(Real r, Real z) const;
    Real Lambda_z(Real r, Real z) const;

};

#endif //GHZ_NUMERIC_KERRMETRIC_HPP
