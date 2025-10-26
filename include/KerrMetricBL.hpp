//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_KERRMETRICBL_HPP
#define GHZ_NUMERIC_KERRMETRICBL_HPP

#pragma once
#include "Metric.hpp"
#include "KerrMetric.hpp"
#include "KerrParams.hpp"

class KerrMetricBL : public Metric {
private:
    KerrParams params;
    KerrMetric kerr_metric;

public:
    explicit KerrMetricBL(const KerrParams& p, const KerrMetric& km);

    //double Sigma(Real r, Real th) const;
    //double Delta(Real r) const;
    //double Lambda(Real r) const;

    [[nodiscard]] std::array<double, 10> g(Real t, Real r, Real th, Real phi) const ;
    [[nodiscard]] std::array<double, 10> ginv(Real t, Real r, Real th, Real phi) const ;
};

#endif //  //GHZ_NUMERIC_KERRMETRICBL_HPP