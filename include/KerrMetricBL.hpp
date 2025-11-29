//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRICBL_HPP
#define GHZ_NUMERIC_KERRMETRICBL_HPP

#pragma once

#include "Metric.hpp"
#include "KerrMetric.hpp"
#include "Coords.hpp"
#include "KerrParams.hpp"

#include <optional>


/** \class KerrMetricBL
* \brief Boyer–Lindquist coordinate representation of the Kerr metric.
*
* This class provides a BL–chart view of a Kerr spacetime.
* It stores the physical Kerr parameters (`M`, `a`) and offers
* metric components `g`, inverse metric `ginv`, and the standard
* Newman–Penrose quantity `rho(r,z)` evaluated in BL coordinates.
*

* The `build()` function precomputes BL-dependent quantities at a
* given coordinate point, and `g()`, `ginv()`, and `rho()` provide
* access to the metric data evaluated at that point.
*/
class KerrMetricBL : public Metric {
private:
    KerrParams params; // parameters of the Kerr metric
    KerrMetric kerr_metric; // reference Kerr metric with all derived parameters computed
    // Cached values computed by build()
    Real sig_, del_, s1_, s2_, a_, M_;

    // Track which point was used for build()
    // in header:

    std::optional<BLCoords> X_cached_;
    bool cache_valid_ = false;

public:
    explicit KerrMetricBL(const KerrParams& p, const KerrMetric& km);

    /// Precompute all BL-dependent quantities at a point.
    /// Must be called before g(), ginv(), rho() with the same coordinates.
    void build(const BLCoords Xbl);

    [[nodiscard]] ghz::SymmetricMatrix4 g(const BLCoords Xbl) const; // metric components g_{\mu\nu}
    [[nodiscard]] ghz::SymmetricMatrix4 ginv(const BLCoords Xbl) const; // inverse metric components g^{\mu\nu}
    [[nodiscard]] ghz::Complex rho(const Real &r, const Real& z) const; // Newman–Penrose quantity rho(r,z) in BL coords
};

#endif //  //GHZ_NUMERIC_KERRMETRICBL_HPP