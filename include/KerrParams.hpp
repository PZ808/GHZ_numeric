//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_KERRPARAMS_HPP
#define GHZ_NUMERIC_KERRPARAMS_HPP

#pragma once
#include <cassert>
#include "GhzTypes.hpp"

using namespace teuk;

struct KerrParams {
    Real M;
    Real a;

    KerrParams(Real mass, Real spin);
};

#endif //GHZ_NUMERIC_KERRPARAMS_HPP
