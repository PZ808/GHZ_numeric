//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRICINGOING_HPP
#define GHZ_NUMERIC_KERRMETRICINGOING_HPP

#pragma once
#include "Metric.hpp"
#include "KerrMetric.hpp"
#include "KerrParams.hpp"

class KerrMetricIngoing : public Metric {
private:
    KerrParams params;
    KerrMetric kerr_metric;

public:
    explicit KerrMetricIngoing(const KerrParams& p, const KerrMetric& km);

    [[nodiscard]] std::array<double, 10> g(double v, double r, double th, double phi_in) const ;
    [[nodiscard]] std::array<double, 10> ginv(double v, double r, double th, double phi_in) const ;
};

#endif //GHZ_NUMERIC_KERRMETRICINGOING_HPP
