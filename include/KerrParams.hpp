//
// Created by Peter Zimmerman on 24.10.25.
//

#ifndef GHZ_NUMERIC_KERRPARAMS_HPP
#define GHZ_NUMERIC_KERRPARAMS_HPP

#pragma once
#include <cassert>


struct KerrParams {
    double M;
    double a;

    KerrParams(double mass, double spin);
};

#endif //GHZ_NUMERIC_KERRPARAMS_HPP
