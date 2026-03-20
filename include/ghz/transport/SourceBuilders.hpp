//
// Created by Peter Zimmerman on 05.03.26.
//

#ifndef HS1DATA_DATA_SOURCEBUILDERS_HPP
#define HS1DATA_DATA_SOURCEBUILDERS_HPP
#pragma once
#include "ghz/core/GhzTypes.hpp"
#include "ghz/transport/RK/ODE.hpp"
#include "ghz/core/MathMacros.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/geom/KerrMetricOutgoing.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/transport/RK/Operators.hpp"      // XDerivs, ModesMK, GHPType

namespace ghz::source {

// x_mmbar source: just returns T_ll if present, else zeros
    spectral::SpectralGHPVectorized
    source_xmmb(size_t Nr, size_t Nz, transport::ModesMK modes, ghp::GHPType out_type,
                const spectral::SpectralGHPVectorized *Tll);


// x_nm source: N[x_mmbar] + T_lm
    spectral::SpectralGHPVectorized
    source_xnm(const transport::XDerivs& xmmbar,               // needs X, P_X, EH_X, EH_P_X, etc (whatever your formula uses)
               const KerrMetricOutgoing& metric,
               const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords>& held,
               const ghp::GHPBackgroundFieldsVectorized& ghp,
               ghp::GHPType out_type,
               const spectral::SpectralGHPVectorized* Tlm);

// x_nn source: T_ln + Re(U[x_mmbar]) + Re(V[x_nm])
    spectral::SpectralGHPVectorized
    source_xnn(const transport::XDerivs& xmmbar,
               const transport::XDerivs& xnm,
               const KerrMetricOutgoing& metric,
               const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords>& held,
               const ghp::GHPBackgroundFieldsVectorized& ghp,
               ghp::GHPType out_type,
               const spectral::SpectralGHPVectorized* Tln);


}

#endif //HS1DATA_DATA_SOURCEBUILDERS_HPP
