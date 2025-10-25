//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_METRIC_HPP
#define GHZ_NUMERIC_METRIC_HPP

#pragma once
#include <array>

class Metric {
public:
    virtual ~Metric() = default;
    //std::array<double, 10> g(double time, double radial, double polar, double azi) const = 0;
    //std::array<double, 10> ginv(double time, double radial, double polar, double azi) const = 0;
};

#endif //GHZ_NUMERIC_METRIC_HPP
