//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
#define GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
#pragma once

#include "Tetrads.hpp"

// family of Kinnersly tetrad classes — one for every possible CoordT type
template <typename CoordT>
class KinnersleyTetrad : public Tetrad {

private:
    std::optional<CoordT> X_cached_;
    bool cache_valid_ = false;

public:

    using Tetrad::Tetrad;
    using Real = teuk::Real;  // or from metric/coords if templated further

    [[nodiscard]] inline Real a()  const { return metric.a(); } // convenience accessors
    [[nodiscard]] inline Real M()  const { return metric.M(); } // convenience accessors
    [[nodiscard]] inline const KerrMetric & get_metric() const { return metric; }

    // Coordinate-dependent tetrad builder
    void build_tetrad_at(const CoordT& X); // builds tetrad at the point X

    SpinCoefficientsGHP get_spin_coeffs_at(const CoordT& X) const;
    WeylScalars         get_weyl_scalars_at(const CoordT& X) const;
    HeldCoefficients get_held_scalars_at(const CoordT& X) const;
    Scalars get_scalars_at(const CoordT& X) const;
    SpinCoefficientsGHP get_ghp_scalars_at(const OutgoingCoordsCompact &X) const;

    Complex rho_at(const Real& r, const Real& z) const;
};


#endif //GHZ_NUMERIC_KINNERSLEYTETRAD_HPP
