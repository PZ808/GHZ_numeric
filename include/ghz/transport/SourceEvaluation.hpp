//
// Created by Peter Zimmerman on 14.03.26.
//

#ifndef HS1DATA_DATA_SOURCEEVALUATION_HPP
#define HS1DATA_DATA_SOURCEEVALUATION_HPP
#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/transport/ODE.hpp"

#include <vector>

namespace ghz::source_eval {

    using GHPSpectral = spectral::SpectralGHPVectorized;
    using RVector = std::vector<teuk::Real>;
    using Real = teuk::Real;
    using Complex = teuk::Complex;
    using StateVec = std::vector<Real>;

    Complex eval_field_slice_linear(const GHPSpectral& field,
                                    const RVector& r_grid,
                                    size_t iz,
                                    Real r);

    ode::SourceFunc make_source_func_xmmbar(const GHPSpectral& source_field,
                                            const RVector& r_grid,
                                            size_t iz);

    ode::SourceFunc make_source_func_xnm(const GHPSpectral& source_field,
                                         const RVector& r_grid,
                                         size_t iz);

    ode::SourceFunc make_source_func_xnn(const GHPSpectral& source_field,
                                         const RVector& r_grid,
                                         size_t iz);

}

#endif //HS1DATA_DATA_SOURCEEVALUATION_HPP
