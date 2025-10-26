//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRICOUTGOING_HPP
#define GHZ_NUMERIC_KERRMETRICOUTGOING_HPP


#pragma once
#include "Metric.hpp"
#include "KerrMetric.hpp"
#include "KerrParams.hpp"

class KerrMetricOutgoing : public Metric {
private:
    KerrParams params;
    KerrMetric kerr_metric;

public:
    explicit KerrMetricOutgoing(const KerrParams& p, const KerrMetric& km);

    [[nodiscard]] std::array<double, 10> g(Real u, Real r, Real th, Real phi_out) const ;
    [[nodiscard]] std::array<double, 10> ginv(Real u, Real r, Real th, Real phi_out) const ;
};

#endif //GHZ_NUMERIC_KERRMETRICOUTGOING_HPP
