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
    Real lambda_, rho0_, Om_;
    Real sig_, del_, s1_, s2_, a_, M_;

public:
    explicit KerrMetricOutgoing(const KerrParams& p, const KerrMetric& km);
    void build(const OutgoingCoords Xout);

    ghz::SymmetricMatrix4 g(const OutgoingCoords Xout) const;
    ghz::SymmetricMatrix4 ginv(const OutgoingCoords Xout) const;

    void compactify(const OutgoingCoords Xout);

    void build_compact(const OutgoingCoords Xout);

    ghz::SymmetricMatrix4 g_tilde(const OutgoingCoords Xout) const;

    void build_compact_from_outgoing(const OutgoingCoords Xout);

    void build_compact(const OutgoingCoordsCompact Xout);

    ghz::SymmetricMatrix4 g_tilde(const OutgoingCoordsCompact Xout) const;
};

#endif //GHZ_NUMERIC_KERRMETRICOUTGOING_HPP
