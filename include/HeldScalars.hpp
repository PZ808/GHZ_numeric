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
#include "../sand/SpectralGHPField.hpp"
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
    HeldField rhopH, rhopH_bar, tauH, tauH_bar, PsiH, PsiH_bar, OmH, OmH_bar;
    HeldBackgroundFields(int Nz) :
            rhopH(Nz, HeldScalar(teuk::zeroC,-2,-2), -2,-2),
            rhopH_bar(Nz, HeldScalar(teuk::zeroC,-2,-2), -2,-2),
            tauH (Nz, HeldScalar(teuk::zeroC,-1,-3), -1,-3),
            tauH_bar(Nz, HeldScalar(teuk::zeroC,-3,-1), -3,-1),
            PsiH (Nz, HeldScalar(teuk::zeroC,-3,-3), -3,-3),
            PsiH_bar(Nz, HeldScalar(teuk::zeroC,-3,-3), -3,-3),
            OmH  (Nz, HeldScalar(teuk::zeroC,-1,-1), -1,-1),
            OmH_bar(Nz,HeldScalar(teuk::zeroC,-1,-1), -1,-1) {}
};

template <typename TetradType, typename CoordType>
HeldBackgroundFields build_held_fields_at(TetradType& tetrad,
                                       const std::vector<teuk::Real>& z_nodes,
                                       CoordType &X)
{
    const int Nz = static_cast<int>(z_nodes.size());
    HeldBackgroundFields held_fields(Nz);

    for (int iz=0; iz<Nz; ++iz) {
        // set z to current collocation point, r can be fixed (e.g. r=const slice)
        X.x2 = z_nodes[iz];

        tetrad.build_tetrad(X);   // build tetrad at this point

        auto scalars = tetrad.get_scalars_at(X);  // get all scalars at the current point X = (r,z)
        // assign directly to HeldFields
        held_fields.rhopH(iz)     = scalars.held_scalars.rhopH;
        held_fields.rhopH_bar(iz) = scalars.held_scalars.rhopH;
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
 * // One can also perform element-wise operations with other HeldFieldVectorized
 * held_fields.PsiH += held_fields.PsiH; // doubles all PsiH values
 * \endcode
 *
 * \note OpenMP is used inside build_held_fields_vectorized to fill all z-points
 *       in parallel. Make sure OpenMP is enabled in your compiler.
 */
namespace ghp {
    struct HeldBackgroundFieldsVectorized {
        HeldFieldVectorized rhopH, rhopH_bar, tauH, tauH_bar, PsiH, PsiH_bar, OmH, OmH_bar;
        // bundle of 1D arrays (along z) for the background Held/GHP scalars
        explicit HeldBackgroundFieldsVectorized(int Nz)
                : rhopH(Nz, HeldScalar(teuk::zeroC, -2, -2), -2, -2),
                  rhopH_bar(Nz, HeldScalar(teuk::zeroC, -2, -2), -2, -2),
                  tauH(Nz, HeldScalar(teuk::zeroC, -1, -3), -1, -3),
                  tauH_bar(Nz, HeldScalar(teuk::zeroC, -3, -1), -3, -1),
                  PsiH(Nz, HeldScalar(teuk::zeroC, -3, -3), -3, -3),
                  PsiH_bar(Nz, HeldScalar(teuk::zeroC, -3, -3), -3, -3),
                  OmH(Nz, HeldScalar(teuk::zeroC, -1, -1), -1, -1),
                  OmH_bar(Nz, HeldScalar(teuk::zeroC, -1, -1), -1, -1)
                  {}
    };

    struct HeldRicciEquationsTypeDVectorized {
        HeldFieldVectorized EdthH_tauH, EdthH_tauH_bar, EdthH_PsiH, EdthH_PsiH_bar,
                EdthH_OmH, EdthH_OmH_bar, EdthH_rhopH;
        HeldFieldVectorized EdthHbar_tauH, EdthHbar_tauH_bar, EdthHbar_PsiH, EdthHbar_PsiH_bar,
                EdthHbar_OmH, EdthHbar_OmH_bar, EdthHbar_rhopH;
        HeldFieldVectorized ThPHeld_rhopH, ThPHeld_tauH, ThPHeld_tauH_bar, ThPHeld_PsiH, ThPHeld_PsiH_bar,
                ThPHeld_OmH, ThPHeld_OmH_bar;
        explicit HeldRicciEquationsTypeDVectorized(int Nz) :
            // Edth = \tilde\edth = {0, -2}
            EdthH_tauH(Nz, GHPScalar(teuk::zeroC, -1, -5), -1, -5),
            EdthH_tauH_bar(Nz, HeldScalar(teuk::zeroC, -3, -3), -3, -3),
            EdthH_PsiH(Nz, HeldScalar(teuk::zeroC, -3, -5), -3, -5),
            EdthH_PsiH_bar(Nz, HeldScalar(teuk::zeroC, -3, -5), -3, -5),
            EdthH_OmH(Nz, HeldScalar(teuk::zeroC, -1, -3), -1, -3),
            EdthH_OmH_bar(Nz, HeldScalar(teuk::zeroC, -1, -3), -1, -3),
            EdthH_rhopH(Nz, HeldScalar(teuk::zeroC, -2, -4), -2, -4),
            // EdthHbar = \tilde\edth' = {-2,0}
            EdthHbar_tauH(Nz, HeldScalar(teuk::zeroC, -3, -3), -3, -3),
            EdthHbar_tauH_bar(Nz, HeldScalar(teuk::zeroC, -5, -1), -5, -1),
            EdthHbar_PsiH(Nz, HeldScalar(teuk::zeroC, -5, -3), -5, -3),
            EdthHbar_PsiH_bar(Nz, HeldScalar(teuk::zeroC, -5, -3), -5, -3),
            EdthHbar_OmH(Nz, HeldScalar(teuk::zeroC, -3, -1), -3, -1),
            EdthHbar_OmH_bar(Nz, HeldScalar(teuk::zeroC, -3, -1), -3, -1),
            EdthHbar_rhopH(Nz, HeldScalar(teuk::zeroC, -4, -2), -4, -2),
            // \tilde\thorn'  = {-1,-1}
            ThPHeld_tauH(Nz, HeldScalar(teuk::zeroC, -2, -4), -2, -4),
            ThPHeld_tauH_bar(Nz, HeldScalar(teuk::zeroC, -4, -2), -4, -2),
            ThPHeld_PsiH(Nz, HeldScalar(teuk::zeroC, -4, -4), -4, -4),
            ThPHeld_PsiH_bar(Nz, HeldScalar(teuk::zeroC, -4, -4), -4, -4),
            ThPHeld_OmH(Nz, HeldScalar(teuk::zeroC, -2, -2), -2, -2),
            ThPHeld_OmH_bar(Nz, HeldScalar(teuk::zeroC, -2, -2), -2, -2),
            ThPHeld_rhopH(Nz, HeldScalar(teuk::zeroC, -3, -3), -3, -3)
        {}
    };

    /// Bundle for background Held fields + their Ricci-type-D combinations
    struct HeldBackgroundAndRicciVectorized {
        HeldBackgroundFieldsVectorized fields;
        HeldRicciEquationsTypeDVectorized riccis;

        explicit HeldBackgroundAndRicciVectorized(int Nz)
                : fields(Nz), riccis(Nz) {}
    };

    // build the Held background fields in vectorized form for all z-nodes in parallel using OpenMP
    template<typename TetradType, typename CoordType>
    HeldBackgroundAndRicciVectorized build_held_fields_vectorized(TetradType &tetrad,
                                                                const std::span<const teuk::Real> z_nodes,
                                                                CoordType &X)
    {
        const int Nz = static_cast<int>(z_nodes.size());
        HeldBackgroundAndRicciVectorized result(Nz);


        auto &held_fields = result.fields;
        auto &held_riccis = result.riccis;

        // Use OpenMP to parallelize over z-nodes
#pragma omp parallel for default(none) firstprivate(tetrad, X, Nz) shared(z_nodes, held_fields, held_riccis)
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

            // create an auto functor for a Held scalar zero with correct weights
                auto zero_held_scalar = [&](int p, int q) {
                    return HeldScalar(Complex(0.0, 0.0), p, q);
                };

            Complex half = Complex(1.0/2.0,0);
            // Type D Ricci equations in Held form using the background fields
            // as in Held Comm. math. Phys. 37, 311—326 (1974) Eqs 5.12-5.21
            // these should hold in vacuum type D spacetime and can be
            // used to verify the correctness of the background fields
            /*
             * Edth on Held scalars from Ricci type-D
             */
            held_riccis.EdthH_tauH(iz) = zero_held_scalar(-1, -5);
            held_riccis.EdthH_tauH_bar(iz) = -half*held_fields.OmH(iz)*(held_fields.rhopH(iz)
                                                                        +held_fields.rhopH_bar(iz))
                                             - half*(held_fields.PsiH(iz)-held_fields.PsiH_bar(iz));
            held_riccis.EdthH_PsiH(iz) = zero_held_scalar(-3, -5);
            held_riccis.EdthH_PsiH_bar(iz) = zero_held_scalar(-3, -5);
            held_riccis.EdthH_OmH(iz) = Complex(2.0,0)*held_fields.tauH(iz);
            held_riccis.EdthH_OmH_bar(iz) = -Complex(2.0,0)*held_fields.tauH(iz);
            held_riccis.EdthH_rhopH(iz) = zero_held_scalar(-2, -4);
            /*
             * EdthBar on Held scalars from Ricci type-D
             */
            held_riccis.EdthHbar_tauH(iz) = -half*held_fields.OmH_bar(iz)*(held_fields.rhopH(iz)
                                                                           +held_fields.rhopH_bar(iz))
                                              - half*(held_fields.PsiH(iz)-held_fields.PsiH_bar(iz) );
            held_riccis.EdthHbar_tauH_bar(iz) = zero_held_scalar(-5, -1);
            held_riccis.EdthHbar_PsiH(iz) = zero_held_scalar(-5, -3);
            held_riccis.EdthHbar_PsiH_bar(iz) = zero_held_scalar(-5, -3);
            held_riccis.EdthHbar_OmH(iz) = -Complex(2.0,0)*held_fields.tauH_bar(iz);
            held_riccis.EdthHbar_OmH_bar(iz) = Complex(2.0,0)*held_fields.tauH_bar(iz);
            held_riccis.EdthHbar_rhopH(iz) = zero_held_scalar(-4, -2);
            /*
             * ThPHeld on Held scalars from Ricci type-D
             */
            held_riccis.ThPHeld_tauH(iz) = zero_held_scalar(-2, -4);
            held_riccis.ThPHeld_tauH_bar(iz) = zero_held_scalar(-4,-2);
            held_riccis.ThPHeld_PsiH(iz) = zero_held_scalar(-4, -4);
            held_riccis.ThPHeld_PsiH_bar(iz) = zero_held_scalar(-4, -4);
            held_riccis.ThPHeld_OmH(iz) = zero_held_scalar(-2, -2);
            held_riccis.ThPHeld_OmH_bar(iz) = zero_held_scalar(-2, -2);
            held_riccis.ThPHeld_rhopH(iz) = zero_held_scalar(-3, -3);
        }

        return result;
    }
} // ghp

#endif //GHZ_NUMERIC_HELDSCALARS_HPP
