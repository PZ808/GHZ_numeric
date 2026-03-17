//
// Created by Peter Zimmerman on 08.03.26.
//

#ifndef HS1DATA_DATA_DOMAIN_HPP
#define HS1DATA_DATA_DOMAIN_HPP

#include "ghz/core/GHZTypes.hpp"


namespace ghz::numeric {

    struct Domain {
        teuk::Real r_lower;
        teuk::Real r_upper;

        Domain(teuk::Real r_lower, teuk::Real r_upper)
                : r_lower(r_lower), r_upper(r_upper) {}
    };
}

#endif //HS1DATA_DATA_DOMAIN_HPP
