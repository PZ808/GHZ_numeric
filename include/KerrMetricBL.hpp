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



class KerrMetricBL : public Metric {
private:
    KerrParams params;
    KerrMetric kerr_metric;
    Real sig_, del_, s1_, s2_, a_, M_;

public:
    explicit KerrMetricBL(const KerrParams& p, const KerrMetric& km);
    void build(const BLCoords Xbl);

    //Real Sigma(Real r, Real th) const;
    //Real Delta(Real r) const;
    //Real Lambda(Real r) const;

    [[nodiscard]] ghz::SymmetricMatrix4 g(const BLCoords Xbl) const;
    [[nodiscard]] ghz::SymmetricMatrix4 ginv(const BLCoords Xbl) const;
};

#endif //  //GHZ_NUMERIC_KERRMETRICBL_HPP