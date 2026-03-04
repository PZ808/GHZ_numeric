//
// Created by Peter Zimmerman on 04.03.26.
//
#ifndef HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP
#define HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP

#pragma once
#include "GhzTypes.hpp"
#include "SpectralGHPFieldVectorized.hpp"
#include "GHPScalars.hpp"
#include "HeldScalars.hpp"
#include "SpectralDiffer.hpp"
#include "KinnersleyTetrad.hpp"

using namespace  spectral;
using namespace ghp;

template <typename CoordType>
class KinnersleyHeldOperators {
public:
    KinnersleyHeldOperators(const spectral::SpectralDiffer& diff,
                            const ghp::HeldBackgroundFieldsVectorized<CoordType>& bg_helds,
                            const KinnersleyTetrad<CoordType>& kin_tetrad)
            : diff_(diff), bg_helds_(bg_helds), kin_tetrad_(kin_tetrad) {}


    void thornPH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                                 SpectralGHPVectorized::RSlice& out) const;

    void edthH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                               SpectralGHPVectorized::RSlice& out) const;
    void edthBarH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                               SpectralGHPVectorized::RSlice& out) const;

    void edthH_dmat_RSlice_inplace(const SpectralGHPVectorized::RSlice &in_RSlice,
                                   SpectralGHPVectorized::RSlice &out_RSlice) const;

    void edthH_bary_inplace(const SpectralGHPVectorized::RSlice &f,
                            SpectralGHPVectorized::RSlice &out) const;

private:
    const spectral::SpectralDiffer& diff_;
    const HeldBackgroundFieldsVectorized<CoordType>& bg_helds_;
    const KinnersleyTetrad<CoordType>& kin_tetrad_;
    const Real a_ = kin_tetrad_.a();
};

#endif //HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP
