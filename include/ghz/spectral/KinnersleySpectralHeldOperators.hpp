//
// Created by Peter Zimmerman on 04.03.26.
//
#ifndef HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP
#define HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP

#pragma once
#include "ghz/core/GhzTypes.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/geom/KinnersleyTetrad.hpp"

using namespace  spectral;
using namespace ghp;

namespace detail {

    static inline Complex quad_extrapolate_endpoint(
            Real z0,
            Real z1, Complex f1,
            Real z2, Complex f2,
            Real z3, Complex f3)
    {
        // Quadratic Lagrange extrapolation to z0 from points (z1,z2,z3).
        const Real d12 = (z1 - z2);
        const Real d13 = (z1 - z3);
        const Real d21 = (z2 - z1);
        const Real d23 = (z2 - z3);
        const Real d31 = (z3 - z1);
        const Real d32 = (z3 - z2);

        Complex L1 = f1 * ((z0 - z2) * (z0 - z3)) / (d12 * d13);
        Complex L2 = f2 * ((z0 - z1) * (z0 - z3)) / (d21 * d23);
        Complex L3 = f3 * ((z0 - z1) * (z0 - z2)) / (d31 * d32);
        return L1 + L2 + L3;
    }

} // namespace detail

enum class EthKind { Eth, EthBar };

template <typename CoordType>
class KinnersleyHeldOperators {
public:
    KinnersleyHeldOperators(const spectral::SpectralDiffer& diff,
                            const ghp::HeldBackgroundFieldsVectorized<CoordType>& bg_helds,
                            const KinnersleyTetrad<CoordType>& kin_tetrad)
            : diff_(diff), bg_helds_(bg_helds), kin_tetrad_(kin_tetrad) {}

    void edth_core_RSlice_with_extrapolation(const spectral::SpectralGHPVectorized::RSlice& in,
                                             const spectral::SpectralGHPVectorized::RSlice& df_dz,
                                             spectral::SpectralGHPVectorized::RSlice& out,
                                             EthKind kind) const;

    void thornPH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                                 SpectralGHPVectorized::RSlice& out) const;

    void edthH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                               SpectralGHPVectorized::RSlice& out) const;
    void edthBarH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                                  SpectralGHPVectorized::RSlice& out) const;

    void edthH_dmat_inplace_RSliceV(const SpectralGHPVectorized::RSlice &in_RSlice,
                                   SpectralGHPVectorized::RSlice &out_RSlice) const;

    void edthH_bary_inplace_RSliceV(const SpectralGHPVectorized::RSlice &f,
                            SpectralGHPVectorized::RSlice &out) const;

private:
    const spectral::SpectralDiffer& diff_;
    const HeldBackgroundFieldsVectorized<CoordType>& bg_helds_;
    const KinnersleyTetrad<CoordType>& kin_tetrad_;
    const Real a_ = kin_tetrad_.a();

    void
    thornPHr_inplace_RSliceV(const KinnersleyTetrad<OutgoingCoords> &ktet,
                             const SpectralGHPVectorized::RSlice &in_RSlice,
                             const SpectralGHPVectorized::RSlice &dr_in_RSlice,
                             SpectralGHPVectorized::RSlice &out_RSlice) const;
};

#endif //HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP
