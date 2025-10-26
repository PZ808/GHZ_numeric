//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRICOUTGOING_HPP
#define GHZ_NUMERIC_KERRMETRICOUTGOING_HPP


#pragma once
#include "Metric.hpp"
#include "Coords.hpp"
#include "KerrMetric.hpp"
#include "KerrParams.hpp"

class KerrMetricOutgoing : public Metric {
private:
    KerrParams params;
    KerrMetric kerr_metric;
    Real sig_, del_, s1_, s2_, a_, M_;

public:
    explicit KerrMetricOutgoing(const KerrParams& p, const KerrMetric& km);
    void build(const OutgoingCoords Xout);

    //[[nodiscard]] ghz::SymmetricMatrix4 g(Real u, Real r, Real z, Real phi_out) const ;
    //[[nodiscard]] ghz::SymmetricMatrix4 ginv(Real u, Real r, Real z, Real phi_out) const ;

    ghz::SymmetricMatrix4 g(const OutgoingCoords Xout) const;
    ghz::SymmetricMatrix4 ginv(const OutgoingCoords Xout) const;
};

#endif //GHZ_NUMERIC_KERRMETRICOUTGOING_HPP
