//
// Created by Peter Zimmerman on 08.03.26.
//

#ifndef HS1DATA_DATA_SPECTRALCOORDINATEMAPS_HPP
#define HS1DATA_DATA_SPECTRALCOORDINATEMAPS_HPP

// ghz/numeric/CoordinateMap1D.hpp
#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/geom/Domain.hpp"
#include <stdexcept>

namespace ghz::numeric {

    struct AffineMap1D {
        teuk::Real a;
        teuk::Real b;

        AffineMap1D(Domain domain) : a(domain.r_lower), b(domain.r_upper) {
            if (a == b) {
                throw std::runtime_error("AffineMap1D: invalid interval");
            }
        }

        [[nodiscard]] teuk::Real toSpectral(teuk::Real r) const {
            return (teuk::two*r - (a + b)) / (b - a);
        }

        [[nodiscard]] teuk::Real toPhysical(teuk::Real x) const {
            return ((b - a)*(x + teuk::one)) / teuk::two + a;
        }

        [[nodiscard]] teuk::Real dxdr() const {
            return teuk::two / (b - a);
        }

        [[nodiscard]] teuk::Real drdx() const {
            return (b - a) / teuk::two;
        }

        [[nodiscard]] std::vector<teuk::Real>
        toPhysicalNodes(const std::vector<teuk::Real>& x_nodes) const {
            std::vector<teuk::Real> r_nodes(x_nodes.size());
            for (size_t i = 0; i < x_nodes.size(); ++i) {
                r_nodes[i] = toPhysical(x_nodes[i]);
            }
            return r_nodes;
        }
    };
}
#endif //HS1DATA_DATA_SPECTRALCOORDINATEMAPS_HPP