//
// Created by Peter Zimmerman on 08.03.26.
//

#ifndef GHZ_NUMERIC_DATADOMAIN_HPP
#define GHZ_NUMERIC_DATADOMAIN_HPP

#pragma once

#include "ghz/core/GhzTypes.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>

namespace ghz::numeric {

    using Real = teuk::Real;

    struct Domain {
        Real lower{};
        Real upper{};

        Domain() = default;

        Domain(Real lower_in, Real upper_in)
                : lower(lower_in), upper(upper_in) {
            if (!(lower < upper)) {
                throw std::invalid_argument("Domain: require lower < upper");
            }
        }

        [[nodiscard]] Real length() const noexcept { return upper - lower; }

        [[nodiscard]] bool contains(Real x) const noexcept {
            return (lower <= x && x <= upper);
        }

        [[nodiscard]] bool interiorContains(Real x) const noexcept {
            return (lower < x && x < upper);
        }
    };

    enum class RegionType {
        Vacuum,
        Puncture
    };

    struct Subdomain1D {
        Domain interval;
        RegionType type{RegionType::Vacuum};
        std::string name{};

        Subdomain1D() = default;

        Subdomain1D(Domain interval_in,
                    RegionType type_in,
                    std::string name_in)
                : interval(interval_in),
                  type(type_in),
                  name(std::move(name_in)) {}
    };

    /**
     * @brief Radial decomposition adapted to a particle at r = r_p.
     *
     * Splits the full radial interval
     *
     *   [r_inner, r_outer]
     *
     * into
     *
     *   vac_left        = [r_inner, r_min]
     *   puncture_left   = [r_min,   r_p]
     *   puncture_right  = [r_p,     r_max]
     *   vac_right       = [r_max,   r_outer]
     *
     * with strict ordering
     *
     *   r_inner < r_min < r_p < r_max < r_outer.
     */
    struct RadialDomainSplit {
        Domain full;

        Real r_p{};
        Real r_min{};
        Real r_max{};

        Subdomain1D vac_left;
        Subdomain1D puncture_left;
        Subdomain1D puncture_right;
        Subdomain1D vac_right;

        RadialDomainSplit() = default;

        RadialDomainSplit(Real r_inner,
                          Real r_outer,
                          Real r_particle,
                          Real r_min_in,
                          Real r_max_in)
                : full(r_inner, r_outer),
                  r_p(r_particle),
                  r_min(r_min_in),
                  r_max(r_max_in)
        {
            validate_();

            vac_left = Subdomain1D(
                    Domain(r_inner, r_min),
                    RegionType::Vacuum,
                    "vac_left"
            );

            puncture_left = Subdomain1D(
                    Domain(r_min, r_p),
                    RegionType::Puncture,
                    "puncture_left"
            );

            puncture_right = Subdomain1D(
                    Domain(r_p, r_max),
                    RegionType::Puncture,
                    "puncture_right"
            );

            vac_right = Subdomain1D(
                    Domain(r_max, r_outer),
                    RegionType::Vacuum,
                    "vac_right"
            );
        }

        // helpers
        [[nodiscard]] std::vector<Subdomain1D> orderedSubdomains() const {
            return {vac_left, puncture_left, puncture_right, vac_right};
        }

        [[nodiscard]] std::vector<Subdomain1D> punctureSubdomains() const {
            return {puncture_left, puncture_right};
        }

        [[nodiscard]] std::vector<Subdomain1D> vacuumSubdomains() const {
            return {vac_left, vac_right};
        }

        [[nodiscard]] bool isPunctureInterface(Real r) const noexcept {
            return r == r_p;
        }
    private:
        void validate_() const {
            if (!(full.lower < r_min)) {
                throw std::invalid_argument(
                        "RadialDomainSplit: require r_inner < r_min");
            }
            if (!(r_min < r_p)) {
                throw std::invalid_argument(
                        "RadialDomainSplit: require r_min < r_p");
            }
            if (!(r_p < r_max)) {
                throw std::invalid_argument(
                        "RadialDomainSplit: require r_p < r_max");
            }
            if (!(r_max < full.upper)) {
                throw std::invalid_argument(
                        "RadialDomainSplit: require r_max < r_outer");
            }
        }
    };

    struct PunctureTwoDomainSplit {
        Real r_min{};
        Real r_p{};
        Real r_max{};

        Subdomain1D left;
        Subdomain1D right;

        PunctureTwoDomainSplit(Real r_min_in, Real r_particle, Real r_max_in)
                : r_min(r_min_in), r_p(r_particle), r_max(r_max_in)
        {
            if (!(r_min < r_p && r_p < r_max)) {
                throw std::invalid_argument(
                        "PunctureTwoDomainSplit: require r_min < r_p < r_max");
            }

            left = Subdomain1D(
                    Domain(r_min, r_p),
                    RegionType::Puncture,
                    "puncture_left"
            );

            right = Subdomain1D(
                    Domain(r_p, r_max),
                    RegionType::Puncture,
                    "puncture_right"
            );
        }
    };
} // namespace ghz::numeric

#endif //GHZ_NUMERIC_DATADOMAIN_HPP
