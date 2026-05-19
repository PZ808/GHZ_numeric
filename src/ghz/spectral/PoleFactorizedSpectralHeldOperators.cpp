//
// Created by Peter Zimmerman on 10.04.26.
//
//
// PoleFactorizedHeldOperators.cpp
//

#include "ghz/spectral/KinnersleySpectralHeldOperators.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"

#include <cassert>
#include <cmath>
#include <vector>
#include <span>

namespace {

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

    inline int spin_from_pq(const int p, const int q)
    {
        // Assumes integer spin, i.e. p-q is even.
        assert(((p - q) % 2) == 0);
        return (p - q) / 2;
    }

// Replace this with your actual metadata accessor for m.
// For example, this might be `sl.modes_.m`, `sl.modes_.m()`, or `sl.mode_m()`.
    template <typename SliceLike>
    inline int slice_mode_m(const SliceLike& sl)
    {
        return sl.modes_.m; // <-- adjust if needed
    }

    inline Real reduced_eth_Cplus(const int m, const int s, const Real z)
    {
        const bool left  = (m + s < 0);
        const bool right = (m - s > 0);

        Real c = 1.0_r;
        if (left)  c *= (1.0_r - z);
        if (right) c *= (1.0_r + z);
        return c;
    }

    inline Real reduced_eth_Bplus(const int m, const int s, const Real z)
    {
        const bool left  = (m + s < 0);
        const bool right = (m - s > 0);

        Real b = 0.0_r;

        if (right) {
            b += Real(m - s) * (left ? (1.0_r - z) : 1.0_r);
        }
        if (left) {
            b += Real(m + s) * (right ? (1.0_r + z) : 1.0_r);
        }

        return b;
    }

    inline Real reduced_ethbar_Cminus(const int m, const int s, const Real z)
    {
        const bool left  = (m + s > 0);
        const bool right = (m - s < 0);

        Real c = 1.0_r;
        if (left)  c *= (1.0_r - z);
        if (right) c *= (1.0_r + z);
        return c;
    }

    inline Real reduced_ethbar_Bminus(const int m, const int s, const Real z)
    {
        const bool left  = (m + s > 0);
        const bool right = (m - s < 0);

        Real b = 0.0_r;

        if (left) {
            b += Real(m + s) * (right ? (1.0_r + z) : 1.0_r);
        }
        if (right) {
            b += Real(m - s) * (left ? (1.0_r - z) : 1.0_r);
        }

        return b;
    }

    template <typename InRSliceLike, typename DzRSliceLike, typename DiffLike>
    void apply_edthH_red_core(
            const InRSliceLike& in_RSlice,
            const DzRSliceLike& dz_RSlice,
            SpectralGHPVectorized::RSlice& out_RSlice,
            const DiffLike& diff,
            const Real a_spin)
    {
        assert(in_RSlice.size() == out_RSlice.size());
        assert(dz_RSlice.size() == out_RSlice.size());

        const size_t Nz = in_RSlice.size();
        if (Nz == 0) return;

        const int p = in_RSlice[0].p();
        const int q = in_RSlice[0].q();
        const int s = spin_from_pq(p, q);
        const int m = slice_mode_m(in_RSlice);

        const Real omega  = in_RSlice.has_omega_mk() ? in_RSlice.omega_mk() : Real(0);
        const Real aomega = a_spin * omega;

        const Real inv_sqrt2 = 1.0_r / std::sqrt(2.0_r);

        for (size_t i = 0; i < Nz; ++i) {
            const Real z = diff.lgl_nodes()[i];

            const Real Cplus = reduced_eth_Cplus(m, s, z);
            const Real Bplus = reduced_eth_Bplus(m, s, z);

            const Complex diag = Complex(Bplus - aomega * Cplus, 0.0_r);

            out_RSlice[i].value() =
                    inv_sqrt2 * (Cplus * dz_RSlice[i].value() + diag * in_RSlice[i].value());

            // eth : (p,q) -> (p, q-2)
            out_RSlice[i].set_pq(p, q-2);
        }
    }

    template <typename InRSliceLike, typename DzRSliceLike, typename DiffLike>
    void apply_edthBarH_red_core(
            const InRSliceLike& in_RSlice,
            const DzRSliceLike& dz_RSlice,
            SpectralGHPVectorized::RSlice& out_RSlice,
            const DiffLike& diff,
            const Real a_spin)
    {
        assert(in_RSlice.size() == out_RSlice.size());
        assert(dz_RSlice.size() == out_RSlice.size());

        const size_t Nz = in_RSlice.size();
        if (Nz == 0) return;

        const int p = in_RSlice[0].p();
        const int q = in_RSlice[0].q();
        const int s = spin_from_pq(p, q);
        const int m = slice_mode_m(in_RSlice);

        const Real omega  = in_RSlice.has_omega_mk() ? in_RSlice.omega_mk() : Real(0);
        const Real aomega = a_spin * omega;

        const Real inv_sqrt2 = 1.0_r / std::sqrt(2.0_r);

        for (size_t i = 0; i < Nz; ++i) {
            const Real z = diff.lgl_nodes()[i];

            const Real Cminus = reduced_ethbar_Cminus(m, s, z);
            const Real Bminus = reduced_ethbar_Bminus(m, s, z);

            const Complex diag = Complex(-Bminus + aomega * Cminus, 0.0_r);

            out_RSlice[i].value() =
                    inv_sqrt2 * (Cminus * dz_RSlice[i].value() + diag * in_RSlice[i].value());

            // ethbar : (p,q) -> (p-2, q)
            out_RSlice[i].set_pq(p-2, q);
        }
    }

} // end anonymous namespace


// -----------------------------------------------------------------------------
// Reduced edth_H on an RSlice using dz D-matrix
// -----------------------------------------------------------------------------

template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthHRed_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> dz_buf(
            Nz, GHPScalar<Complex>(Complex(0.0_r, 0.0_r), p, q));

    SpectralGHPVectorized::RSlice dz_RSlice(
            dz_buf.data(), Nz, in_RSlice.modes_, in_RSlice.r_index());

    std::span<const GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> dz_span(dz_RSlice.data_ptr, Nz);

    diff_.dz_Dmatrix(in_span, dz_span);

    apply_edthH_red_core(in_RSlice, dz_RSlice, out_RSlice, diff_, a_); // <-- use your spin param member
}

template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthHRed_inplace_RSliceV(
        const SpectralGHPVectorized::ConstRSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> dz_buf(
            Nz, GHPScalar<Complex>(Complex(0.0_r, 0.0_r), p, q));

    SpectralGHPVectorized::RSlice dz_RSlice(
            dz_buf.data(), Nz, in_RSlice.modes_, in_RSlice.r_index());

    std::span<const GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> dz_span(dz_RSlice.data_ptr, Nz);

    diff_.dz_Dmatrix(in_span, dz_span);

    apply_edthH_red_core(in_RSlice, dz_RSlice, out_RSlice, diff_, a_);
}


// -----------------------------------------------------------------------------
// Reduced edthBar_H on an RSlice using dz D-matrix
// -----------------------------------------------------------------------------

template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthBarHRed_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> dz_buf(
            Nz, GHPScalar<Complex>(Complex(0.0_r, 0.0_r), p, q));

    SpectralGHPVectorized::RSlice dz_RSlice(
            dz_buf.data(), Nz, in_RSlice.modes_, in_RSlice.r_index());

    std::span<const GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> dz_span(dz_RSlice.data_ptr, Nz);

    diff_.dz_Dmatrix(in_span, dz_span);

    apply_edthBarH_red_core(in_RSlice, dz_RSlice, out_RSlice, diff_, a_);
}

template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthBarHRed_inplace_RSliceV(
        const SpectralGHPVectorized::ConstRSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> dz_buf(
            Nz, GHPScalar<Complex>(Complex(0.0_r, 0.0_r), p, q));

    SpectralGHPVectorized::RSlice dz_RSlice(
            dz_buf.data(), Nz, in_RSlice.modes_, in_RSlice.r_index());

    std::span<const GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> dz_span(dz_RSlice.data_ptr, Nz);

    diff_.dz_Dmatrix(in_span, dz_span);

    apply_edthBarH_red_core(in_RSlice, dz_RSlice, out_RSlice, diff_, a_);
}


// -----------------------------------------------------------------------------
// Optional barycentric versions
// -----------------------------------------------------------------------------

template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthHRed_bary_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> dz_buf(
            Nz, GHPScalar<Complex>(Complex(0.0_r, 0.0_r), p, q));

    SpectralGHPVectorized::RSlice dz_RSlice(
            dz_buf.data(), Nz, in_RSlice.modes_, in_RSlice.r_index());

    diff_.dz_barycentric_RSlice(in_RSlice, dz_RSlice);

    apply_edthH_red_core(in_RSlice, dz_RSlice, out_RSlice, diff_, a_);
}

template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthBarHRed_bary_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> dz_buf(
            Nz, GHPScalar<Complex>(Complex(0.0_r, 0.0_r), p, q));

    SpectralGHPVectorized::RSlice dz_RSlice(
            dz_buf.data(), Nz, in_RSlice.modes_, in_RSlice.r_index());

    diff_.dz_barycentric_RSlice(in_RSlice, dz_RSlice);

    apply_edthBarH_red_core(in_RSlice, dz_RSlice, out_RSlice, diff_, a_);
}


// -----------------------------------------------------------------------------
// Reduced radial operators = raw radial operators
// -----------------------------------------------------------------------------

template <>
void KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::thornRed_inplace(
        const SpectralGHPVectorized& in,
        SpectralGHPVectorized& out) const
{
    thorn_inplace(in, out);
}

template <>
void KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::thornPHRed_inplace(
        const SpectralGHPVectorized& in,
        SpectralGHPVectorized& out) const
{
    thornPH_inplace(in, out);
}

template <>
void KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::thornPHrRed_inplace(
        const KinnersleyTetrad<OutgoingCoords>& ktet,
        const SpectralGHPVectorized& in,
        SpectralGHPVectorized& out) const
{
    thornPHr_inplace(ktet, in, out);
}


// -----------------------------------------------------------------------------
// Full-field wrappers
// -----------------------------------------------------------------------------

template <>
void KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::edthHRed_inplace(
        const SpectralGHPVectorized& in,
        SpectralGHPVectorized& out) const
{
    assert(in.Nr() == out.Nr() && in.Nz() == out.Nz());

    for (size_t ir = 0; ir < in.Nr(); ++ir) {
        auto in_r  = in.slice_R(ir);
        auto out_r = out.slice_R(ir);
        edthHRed_inplace_RSliceV(in_r, out_r);
    }
}

template <>
void KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::edthBarHRed_inplace(
        const SpectralGHPVectorized& in,
        SpectralGHPVectorized& out) const
{
    assert(in.Nr() == out.Nr() && in.Nz() == out.Nz());

    for (size_t ir = 0; ir < in.Nr(); ++ir) {
        auto in_r  = in.slice_R(ir);
        auto out_r = out.slice_R(ir);
        edthBarHRed_inplace_RSliceV(in_r, out_r);
    }
}