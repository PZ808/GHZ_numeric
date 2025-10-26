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
    [[nodiscard]] ghz::SymmetricMatrix4 g(double v, double r, double th, double ph) const ;
    [[nodiscard]] ghz::SymmetricMatrix4 ginv(double v, double r, double th, double ph) const ;

};

#endif //GHZ_NUMERIC_KERRMETRICINGOING_HPP
