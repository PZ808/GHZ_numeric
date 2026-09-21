//
// AngularProjection.hpp
//
// Projection of sampled m-omega data onto a single angular mode.
//

#ifndef GHZ_SOURCE_ANGULARPROJECTION_HPP
#define GHZ_SOURCE_ANGULARPROJECTION_HPP
#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/source/TeukolskySource.hpp"

#include <vector>

namespace ghz::source {

    using teuk::Real;
    using teuk::Complex;

    struct AngularModeSamples {
        int spin = 0;
        int ell = 0;
        int m = 0;
        Real spheroidicity = Real(0);

        // Samples of _s S_{\ell m}(z; c) at the angular grid nodes and phi=0.
        // The projection uses the complex conjugate of these values.
        std::vector<Complex> values;
    };

    // Computes 2 pi int_{-1}^{1} f_momega(z) conj(S_lm(z; c)) dz from
    // collocation samples and externally supplied quadrature weights.
    [[nodiscard]] Complex project_angular_slice_z(
            const std::vector<Complex>& source_values,
            const AngularModeSamples& harmonic,
            const std::vector<Real>& weights_z);

    // Applies project_angular_slice_z to every radial row of a sampled source.
    [[nodiscard]] std::vector<Complex> project_angular_mode(
            const TeukolskyField& source,
            const AngularModeSamples& harmonic,
            const std::vector<Real>& weights_z);

} // namespace ghz::source

#endif // GHZ_SOURCE_ANGULARPROJECTION_HPP
