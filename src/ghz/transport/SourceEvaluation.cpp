//
// Created by Peter Zimmerman on 14.03.26.
//
#include "ghz/transport/SourceEvaluation.hpp"

namespace ghz::source_eval {

    inline ode::SourceFunc make_source_func_xnm(const GHPSpectral &source_field,
                                                const RVector &r_grid,
                                                size_t iz) {
        return [&source_field, &r_grid, iz](Real r, Real /*z*/) -> StateVec {
            const Complex S = eval_field_slice_linear(source_field, r_grid, iz, r);
            return StateVec{0.0, std::real(S), 0.0, std::imag(S)};
        };
    }

    inline ode::SourceFunc make_source_func_xnn(const GHPSpectral &source_field,
                                                const RVector &r_grid,
                                                size_t iz) {
        return [&source_field, &r_grid, iz](Real r, Real /*z*/) -> StateVec {
            const Complex S = eval_field_slice_linear(source_field, r_grid, iz, r);
            return StateVec{std::real(S), std::imag(S)};
        };
    }

    inline ode::SourceFunc make_source_func_xmmbar(const GHPSpectral &source_field,
                                                   const RVector &r_grid,
                                                   size_t iz) {
        return [&source_field, &r_grid, iz](Real r, Real /*z*/) -> StateVec {
            const Complex S = eval_field_slice_linear(source_field, r_grid, iz, r);
            return StateVec{0.0, std::real(S), 0.0, std::imag(S)};
        };
    }
}