//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRICINGOING_HPP
#define GHZ_NUMERIC_KERRMETRICINGOING_HPP

#pragma once
#include "Metric.hpp"
#include "KerrMetric.hpp"
#include "Coords.hpp"
#include "KerrParams.hpp"

class KerrMetricIngoing : public Metric {
private:
    KerrParams params;
    KerrMetric kerr_metric;
    Real sig_, del_, s1_, s2_, a_, M_;

public:
    explicit KerrMetricIngoing(const KerrParams& p, const KerrMetric& km);
    void build(const OutgoingCoords Xout);

    [[nodiscard]] ghz::SymmetricMatrix4 g(const IngoingCoords Xin) const ;
    [[nodiscard]] ghz::SymmetricMatrix4 ginv(const IngoingCoords Xin) const ;

};

#endif //GHZ_NUMERIC_KERRMETRICINGOING_HPP
