//
// Created by Peter Zimmerman on 24.10.25.
//

#include "ghz/geom/KerrParams.hpp"
#include <cmath>
#include <cassert>
#include "ghz/core/MathMacros.hpp"

KerrParams::KerrParams(Real mass, Real spin) : M(mass), a(spin)  {
    assert(M > 0.0 && "Mass must be positive");
    assert(math::Abs(teuk::Real(a)) <= M && "|a| must be <= M for Kerr black hole");
}