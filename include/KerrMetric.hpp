//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRIC_HPP
#define GHZ_NUMERIC_KERRMETRIC_HPP

#pragma once
#include <array>
#include "Metric.hpp"
#include "KerrParams.hpp"


class KerrMetric : public Metric {
private:
    KerrParams params;


public:
    explicit KerrMetric(const KerrParams& p);
    virtual ~KerrMetric() = default;

    //virtual double Delta(double r) const = 0;
    //virtual double Sigma(double r, double th) const = 0;


    double M() const;
    double a() const;

    double r_plus() const;
    double r_minus() const;

    double Om_plus() const;
    double Om_minus() const;

    double kappa_plus() const;
    double kappa_minus() const;


    virtual double Sigma(double r, double theta) const;
    virtual double  Delta(double r) const;
    virtual double Lambda(double r, double theta) const;


    std::array<double, 10> g(double t, double r, double th, double ph) const ;
    std::array<double, 10> ginv(double t, double r, double th, double ph) const ;
};

#endif //GHZ_NUMERIC_KERRMETRIC_HPP
