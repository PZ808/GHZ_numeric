//
// Created by Peter Zimmerman on 25.10.25.
//

//
// Created by Peter Zimmerman on 25.10.25.
//

#include "../include/WeylScalars.hpp"

#include <sstream>
#include <cmath>

using Complex = teuk::Complex;

void WeylScalars::set(WeylScalarType type, Complex value) {
    weyl_scalars[type] = value;
}

Complex WeylScalars::get(WeylScalarType type) const {
    auto it = weyl_scalars.find(type);
    if(it != weyl_scalars.end()) return it->second;
    return teuk::zeroC;
}