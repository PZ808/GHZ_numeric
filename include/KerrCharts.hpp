//
// Created by Peter Zimmerman on 26.11.25.
//

#ifndef GHZ_NUMERIC_KERRCHARTS_HPP
#define GHZ_NUMERIC_KERRCHARTS_HPP

#include "KerrMetric.hpp"
#include "KerrMetricBL.hpp"
#include "KerrMetricIngoing.hpp"
#include "KerrMetricOutgoing.hpp"
#include "Coords.hpp"

// ------------------------------------------------------------
// BLChart
// ------------------------------------------------------------
struct BLChart {
    using Coord = BLCoords;
    using CompactCoord = void; // not used
    using MetricType = KerrMetricBL; // <- key

    static void build_at(KerrMetricBL & km, const Coord& X) { km.build_at(X); }
    static teuk::SymmetricMatrix4 g(KerrMetricBL & km, const Coord& X) { return km.g(X); }
    static teuk::SymmetricMatrix4 ginv(KerrMetricBL& km, const Coord& X) { return km.ginv(X); }
};

// ------------------------------------------------------------
// OutgoingChart
// ------------------------------------------------------------
struct OutgoingChart {
    using Coord = OutgoingCoords;
    using CompactCoord = OutgoingCoordsCompact;

    static void build_at(KerrMetricOutgoing& km, const Coord& X) {
        km.build_at(X);
    }

    static teuk::SymmetricMatrix4 g(KerrMetricOutgoing& km, const Coord& X) { return km.g(X); }
    static teuk::SymmetricMatrix4 ginv(KerrMetricOutgoing& km, const Coord& X) { return km.ginv(X); }
    static teuk::SymmetricMatrix4 g_tilde(const KerrMetricOutgoing& km, const CompactCoord& Xc) { return km.g_tilde(Xc); }
};

// ------------------------------------------------------------
// IngoingChart
// ------------------------------------------------------------
struct IngoingChart {
    using Coord = IngoingCoords;
    using CompactCoord = void;

    static void build_at(KerrMetricIngoing& km, const Coord& X) { km.build_at(X); }
    static teuk::SymmetricMatrix4 g(KerrMetricIngoing& km, const Coord& X) { return km.g(X); }
    static teuk::SymmetricMatrix4 ginv(KerrMetricIngoing& km, const Coord& X) { return km.ginv(X); }
};


#endif //GHZ_NUMERIC_KERRCHARTS_HPP
