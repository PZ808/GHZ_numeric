//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRICBL_HPP
#define GHZ_NUMERIC_KERRMETRICBL_HPP

#pragma once
#include "Metric.hpp"
#include "KerrMetric.hpp"
#include "KerrParams.hpp"

class KerrMetricBL : public Metric {
private:
    KerrParams params;
    KerrMetric kerr_metric;

public:
    explicit KerrMetricBL(const KerrParams& p, const KerrMetric& km);

    //double Sigma(double r, double th) const;
    //double Delta(double r) const;
    //double Lambda(double r) const;

    [[nodiscard]] std::array<double, 10> g(double t, double r, double th, double phi) const ;
    [[nodiscard]] std::array<double, 10> ginv(double t, double r, double th, double phi) const ;
};

#endif //  //GHZ_NUMERIC_KERRMETRICBL_HPP