//
// Created by Peter Zimmerman on 25.10.25.
//

#ifndef GHZ_NUMERIC_HELDSCALARS_HPP
#define GHZ_NUMERIC_HELDSCALARS_HPP

#pragma once
#include "GHPScalars.hpp"
#include "WeylScalars.hpp"
#include "SpinCoeffsNP.hpp"
#include "MathMacros.hpp"
#include "GHPFieldVectorized.hpp"
#include "SpectralGHPField.hpp"
#include <complex>
#include <string>


using Complex = teuk::Complex;


class HeldScalar : public GHPScalar<Complex> {
public:
    HeldScalar() = default;
    HeldScalar(const Complex& v, int p, int q) : GHPScalar<Complex>(v, p, q) {}

};

class HeldField : public GHPField {
public:
    HeldField() = default;
    HeldField(int Nz,
              Scalar init = Scalar(teuk::zeroC, 0, 0), int p = 0, int q = 0)
            :     GHPField(1, Nz, init, p, q) {}

    Scalar& operator()(int iz) {
        return GHPField::operator()(0, iz);
    }
    const Scalar& operator()(int iz) const {
        return GHPField::operator()(0, iz);
    }

};

struct HeldCoefficients {
    HeldCoefficients(const SpinCoefficientsGHP &sc_ghp, const WeylScalars &weyl_scs);

    HeldScalar rhopH, tauH,  rhopH_bar, tauH_bar;
    HeldScalar PsiH, PsiH_bar;
    HeldScalar OmH, OmH_bar;

    HeldCoefficients() = default;

};

struct HeldBackgroundFields {
    HeldField rhopH, tauH, tauH_bar, PsiH, PsiH_bar, OmH, OmH_bar;
    HeldBackgroundFields(int Nz) :
            rhopH(Nz, HeldScalar(teuk::zeroC,-2,-2), -2,-2),
            tauH (Nz, HeldScalar(teuk::zeroC,-1,-3), -1,-3),
            tauH_bar(Nz, HeldScalar(teuk::zeroC,-3,-1), -3,-1),
            PsiH (Nz, HeldScalar(teuk::zeroC,-3,-3), -3,-3),
            PsiH_bar(Nz, HeldScalar(teuk::zeroC,-3,-3), -3,-3),
            OmH  (Nz, HeldScalar(teuk::zeroC,-1,-1), -1,-1),
            OmH_bar(Nz,HeldScalar(teuk::zeroC,-1,-1), -1,-1) {}
};

template <typename TetradType, typename CoordType>
HeldBackgroundFields build_held_fields(TetradType& tetrad,
                                       const std::vector<teuk::Real>& z_nodes,
                                       CoordType &X)
{
    const int Nz = static_cast<int>(z_nodes.size());
    HeldBackgroundFields held_fields(Nz);

    for (int iz = 0; iz < Nz; ++iz) {
        X.x2 = z_nodes[iz];       // set current z to LGL collocation point

        tetrad.build_tetrad(X);   // build tetrad at this point

        auto scalars = tetrad.get_scalars_at(X);  // get all scalars
        // assign directly to HeldFields
        held_fields.rhopH(iz)     = scalars.held_scalars.rhopH;
        held_fields.tauH(iz)      = scalars.held_scalars.tauH;
        held_fields.tauH_bar(iz)  = scalars.held_scalars.tauH_bar;
        held_fields.PsiH(iz)      = scalars.held_scalars.PsiH;
        held_fields.PsiH_bar(iz)  = scalars.held_scalars.PsiH_bar;
        held_fields.OmH(iz)       = scalars.held_scalars.OmH;
        held_fields.OmH_bar(iz)   = scalars.held_scalars.OmH_bar;
    }

    return held_fields;
}
/*!
 * \struct HeldBackgroundFieldsVectorized
 * \brief Holds 1D background GHP fields along the z-direction in vectorized storage
 *
 * Each member is a HeldFieldVectorized containing a single GHP scalar
 * component (rhopH, tauH, tauH_bar, PsiH, PsiH_bar, OmH, OmH_bar).
 * The fields are stored in a contiguous vector for cache-friendly access
 * and support OpenMP parallel operations.
 *
 * Usage example:
 * \code
 * // LGL collocation points along z
 * std::vector<teuk::Real> z_nodes = {...};
 *
 * // Initial coordinates and tetrad
 * CoordType X0;
 * TetradType tetrad;
 *
 * // Build vectorized background fields
 * auto held_fields = ghp::build_held_fields_vectorized(tetrad, z_nodes, X0);
 *
 * // Accessing a scalar at index 3
 * auto psi_val = held_fields.PsiH(3);
 * std::cout << held_fields.PsiH.str(3);
 *
 * // You can also perform element-wise operations with other HeldFieldVectorized
 * held_fields.PsiH += held_fields.PsiH; // doubles all PsiH values
 * \endcode
 *
 * \note OpenMP is used inside build_held_fields_vectorized to fill all z-points
 *       in parallel. Make sure OpenMP is enabled in your compiler.
 */

namespace ghp {
    struct HeldBackgroundFieldsVectorized {
        HeldFieldVectorized rhopH, tauH, tauH_bar, PsiH, PsiH_bar, OmH, OmH_bar;

        explicit HeldBackgroundFieldsVectorized(int Nz)
                : rhopH(Nz, HeldScalar(teuk::zeroC, -2, -2), -2, -2),
                  tauH(Nz, HeldScalar(teuk::zeroC, -1, -3), -1, -3),
                  tauH_bar(Nz, HeldScalar(teuk::zeroC, -3, -1), -3, -1),
                  PsiH(Nz, HeldScalar(teuk::zeroC, -3, -3), -3, -3),
                  PsiH_bar(Nz, HeldScalar(teuk::zeroC, -3, -3), -3, -3),
                  OmH(Nz, HeldScalar(teuk::zeroC, -1, -1), -1, -1),
                  OmH_bar(Nz, HeldScalar(teuk::zeroC, -1, -1), -1, -1) {}
    };

    template<typename TetradType, typename CoordType>
    HeldBackgroundFieldsVectorized build_held_fields_vectorized(TetradType &tetrad,
                                                                const std::span<const teuk::Real> z_nodes,
                                                                CoordType &X)
    {
        const int Nz = static_cast<int>(z_nodes.size());
        HeldBackgroundFieldsVectorized held_fields(Nz);

        // Use OpenMP to parallelize over z-nodes
#pragma omp parallel for default(none) firstprivate(tetrad, X, Nz) shared(z_nodes, held_fields)
        for (int iz = 0; iz < Nz; ++iz) {

            CoordType X_local = X;   // make a local copy for thread safety
            TetradType tetrad_local = tetrad; // copy tetrad for this thread
            X_local.x2 = z_nodes[iz]; // set current z to collocation point
            tetrad_local.build_tetrad_at(X_local);   // build tetrad at this point
            auto scalars = tetrad_local.get_scalars_at(X_local); // get all scalars

            // assign directly to HeldFieldVectorized
            held_fields.rhopH(iz) = scalars.held_scalars.rhopH;
            held_fields.tauH(iz) = scalars.held_scalars.tauH;
            held_fields.tauH_bar(iz) = scalars.held_scalars.tauH_bar;
            held_fields.PsiH(iz) = scalars.held_scalars.PsiH;
            held_fields.PsiH_bar(iz) = scalars.held_scalars.PsiH_bar;
            held_fields.OmH(iz) = scalars.held_scalars.OmH;
            held_fields.OmH_bar(iz) = scalars.held_scalars.OmH_bar;
        }

        return held_fields;
    }
} // ghp

#endif //GHZ_NUMERIC_HELDSCALARS_HPP
