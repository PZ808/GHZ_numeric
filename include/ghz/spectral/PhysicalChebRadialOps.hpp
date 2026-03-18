//
// Created by Peter Zimmerman on 18.03.26.
//

#ifndef HS1DATA_DATA_PHYSICALCHEBRADIALOPS_HPP
#define HS1DATA_DATA_PHYSICALCHEBRADIALOPS_HPP
#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/geom/DataDomain.hpp"
#include "ghz/spectral/SpectralCoordinateMaps.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/ghp/GHPScalars.hpp"

#include <boost/numeric/ublas/operation.hpp>
#include <span>
#include <vector>
#include <stdexcept>

namespace ghz::numeric {

    namespace ublas = boost::numeric::ublas;

    struct PhysicalChebRadialOps {
        using Real    = teuk::Real;
        using Complex = teuk::Complex;
        using matrixR = spectral::matrix<Real>;
        using RVector = spectral::RVector;
        using CVector = spectral::CVector;

        AffineMap1D map;

        RVector x_nodes;   // reference nodes in [-1,1]
        RVector r_nodes;   // physical nodes in [a,b]
        RVector wx;        // barycentric weights on x-nodes

        matrixR Dx;        // d/dx
        matrixR Dxx;       // d^2/dx^2
        matrixR Dr;        // d/dr
        matrixR Drr;       // d^2/dr^2

        explicit PhysicalChebRadialOps(
                const spectral::SpectralDiffer& differ,
                Domain domain)
                : map(domain),
                  x_nodes(differ.cl_nodes()),
                  r_nodes(map.toPhysicalNodes(differ.cl_nodes())),
                  wx(differ.x_weights()),
                  Dx(differ.Dx_matrix()),
                  Dxx(ublas::prod(differ.Dx_matrix(), differ.Dx_matrix())),
                  Dr(map.dxdr() * differ.Dx_matrix()),
                  Drr((map.dxdr() * map.dxdr()) *
                      ublas::prod(differ.Dx_matrix(), differ.Dx_matrix()))
        {}

        [[nodiscard]] Real toSpectral(Real r) const { return map.toSpectral(r); }
        [[nodiscard]] Real toPhysical(Real x) const { return map.toPhysical(x); }

        [[nodiscard]] const RVector& x()   const { return x_nodes; }
        [[nodiscard]] const RVector& r()   const { return r_nodes; }
        [[nodiscard]] const matrixR& Dx_matrix()  const { return Dx; }
        [[nodiscard]] const matrixR& Dxx_matrix() const { return Dxx; }
        [[nodiscard]] const matrixR& Dr_matrix()  const { return Dr; }
        [[nodiscard]] const matrixR& Drr_matrix() const { return Drr; }

        // -------------------------------------------------------------
        // Interpolate a fixed-z slice given as complex nodal values on x-nodes.
        // Returns f(r0) and df/dr(r0).
        // -------------------------------------------------------------
        [[nodiscard]] std::pair<Complex, Complex>
        interp_value_and_dr_from_complex_values(
                const CVector& fvals,
                Real r0) const
        {
            if (fvals.size() != x_nodes.size()) {
                throw std::runtime_error(
                        "interp_value_and_dr_from_complex_values: size mismatch.");
            }

            const Real x0 = map.toSpectral(r0);

            auto out = spectral::SpectralDiffer::barycentric_interp_and_derivative(
                    x_nodes, fvals, wx, x0);

            // barycentric routine returns df/dx; convert to df/dr
            out.second *= map.dxdr();
            return out;
        }

        // -------------------------------------------------------------
        // Same, but input is a GHP slice.
        // -------------------------------------------------------------
        [[nodiscard]] std::pair<GHPScalar<Complex>, GHPScalar<Complex>>
        interp_value_and_dr_from_slice(
                std::span<const GHPScalar<Complex>> f_slice,
                Real r0) const
        {
            if (f_slice.size() != x_nodes.size()) {
                throw std::runtime_error(
                        "interp_value_and_dr_from_slice: size mismatch.");
            }

            CVector vals(f_slice.size());
            for (size_t i = 0; i < f_slice.size(); ++i) {
                vals[i] = f_slice[i].value();
            }

            const auto [v, dvdr] =
                    interp_value_and_dr_from_complex_values(vals, r0);

            return {
                    GHPScalar<Complex>(v,    f_slice[0].p(), f_slice[0].q()),
                    GHPScalar<Complex>(dvdr, f_slice[0].p(), f_slice[0].q())
            };
        }

        // -------------------------------------------------------------
        // Convenience: build complex nodal values for one fixed iz slice.
        // -------------------------------------------------------------
        [[nodiscard]] CVector
        extract_complex_zslice(const spectral::GHPSpectralField& field, size_t iz) const
        {
            const size_t Nr = x_nodes.size();
            if (field.Nr() != Nr) {
                throw std::runtime_error("extract_complex_zslice: Nr mismatch.");
            }
            if (iz >= field.Nz()) {
                throw std::runtime_error("extract_complex_zslice: iz out of range.");
            }

            CVector vals(Nr);
            for (size_t ir = 0; ir < Nr; ++ir) {
                vals[ir] = field.at(ir, iz).value();
            }
            return vals;
        }
    };

} // namespace ghz::numeric
#endif //HS1DATA_DATA_PHYSICALCHEBRADIALOPS_HPP
