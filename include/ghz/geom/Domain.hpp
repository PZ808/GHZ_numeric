//
// Created by Peter Zimmerman on 08.03.26.
//

#ifndef GHZ_NUMERIC_DATA_DOMAIN_HPP
#define GHZ_NUMERIC_DATA_DOMAIN_HPP

#include "ghz/core/GhzTypes.hpp"


namespace ghz::numeric {

    struct RadialDomain {
        teuk::Real r_lower;
        teuk::Real r_upper;


        RadialDomain(teuk::Real r_lower, teuk::Real r_upper)
                : r_lower(r_lower), r_upper(r_upper) {}
    };
}

#endif //GHZ_NUMERIC_DATA_DOMAIN_HPP
