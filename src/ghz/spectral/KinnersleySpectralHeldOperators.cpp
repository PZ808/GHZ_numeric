//
// Created by Peter Zimmerman on 04.03.26.
// KinnersleySpectralHeldOperators.cpp
//

#include "ghz/spectral/KinnersleySpectralHeldOperators.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"

#include <cassert>
#include <cmath>
#include <vector>
#include <span>

/**
 * @name eth_core_RSlice_with_extrapolation
 * @tparam CoordType coordinate type (e.g., Outgoing Boyer-Lindquist, Kerr-Schild)
 * @param in input spectral field on a single RSlice (fixed r, varying z) \n
 * @param df_dz precomputed spectral derivative d/dz of the input field on the same RSlice \n
 * @param out output spectral field where the result will be written, same size as in and df_dz \n
 * @param kind specifies whether to compute the edth or edthBar operator \n
 * @brief Compute the core part of the edth/edthBar operator on a single RSlice, with safe handling of endpoints.
 */
template <>
void KinnersleyHeldOperators<OutgoingCoords>::edth_core_RSlice_with_extrapolation(
        const spectral::SpectralGHPVectorized::RSlice& in,
        const spectral::SpectralGHPVectorized::RSlice& df_dz,
        spectral::SpectralGHPVectorized::RSlice& out,
        EthKind kind) const
{
    assert(in.size() == out.size());
    assert(in.size() == df_dz.size());

    const size_t N = in.size();
    if (N == 0) return;

    // Need at least 4 points for quadratic endpoint extrapolation.
    const bool can_quad = (N >= 4);

    // metadata (assume uniform p,q over slice)
    const int p = in[0].p();
    const int q = in[0].q();
    const Complex s = Complex(Real(p - q) * Real(0.5), 0.0);
    const Real m = Real(in.m());
    const Real aw = a_ * in.omega_mk();

    const int out_p = (kind == EthKind::Eth) ? (p + 0) : (p - 2);
    const int out_q = (kind == EthKind::Eth) ? (q - 2) : (q + 0);

    const Complex pref = Complex(-Real(1.0) / Real(std::sqrt(2.0)), 0.0);

    const auto& z_nodes = diff_.lgl_nodes();
    assert(z_nodes.size() == N);

    // Interior only: i = 1..N-2
#pragma omp parallel for default(none) shared(in, df_dz, out, z_nodes, kind) \
    firstprivate(N, s, m, aw, pref, out_p, out_q)
    for (size_t i = 1; i<N-1; ++i) {
        const Real z = z_nodes[i];
        const Real fac_r = std::sqrt(std::max<Real>(Real(0.0), Real(1.0) - z*z));

        if (fac_r <= Real(0.0)) {
            out[i] = GHPScalar<Complex>(Complex(0.0, 0.0), out_p, out_q);
            continue;
        }

        const Complex factor(fac_r, 0.0);

        const Complex fval = in[i].value();
        const Complex dfz  = df_dz[i].value();

        // Combine the pole-singular terms BEFORE dividing by sin(theta)=factor.
        // m-mode: ∂_φ f = i m f, and your original formulas imply:
        //   Eth:    (+ i*dfphi - s z f)/factor  with dfphi=i m f => i*dfphi = -m f
        //   EthBar: (- i*dfphi + s z f)/factor  with dfphi=i m f => -i*dfphi = +m f
        Complex singular_num;
        if (kind == EthKind::Eth) {
            singular_num = Complex(-m, 0.0) * fval - s * z * fval;
        } else {
            singular_num = Complex(+m, 0.0) * fval + s * z * fval;
        }
        const Complex singular = singular_num / factor;

        // aw term sign differs for Eth vs EthBar in your code
        const Complex aw_term = (kind == EthKind::Eth)
                                ? (Complex(aw, 0.0) * factor * fval)
                                : (Complex(-aw, 0.0) * factor * fval);

        // dz term sign: matches your posted Eth and EthBar implementations.
        // If you later decide EthBar should flip this sign, change here.
        const Complex dz_term = -factor * dfz;

        const Complex val = pref * (dz_term + singular + aw_term);
        out[i] = GHPScalar<Complex>(val, out_p, out_q);
    }

    // Endpoints (0 and N-1): extrapolate from interior
    if (N == 1) {
        out[0] = GHPScalar<Complex>(Complex(0.0, 0.0), out_p, out_q);
        return;
    }

    if (can_quad) {
        out[0].value() = detail::quad_extrapolate_endpoint(
                z_nodes[0],
                z_nodes[1], out[1].value(),
                z_nodes[2], out[2].value(),
                z_nodes[3], out[3].value());
        out[0].set_pq(out_p, out_q);

        out[N - 1].value() = detail::quad_extrapolate_endpoint(
                z_nodes[N - 1],
                z_nodes[N - 2], out[N - 2].value(),
                z_nodes[N - 3], out[N - 3].value(),
                z_nodes[N - 4], out[N - 4].value());
        out[N - 1].set_pq(out_p, out_q);
    } else {
        // N=2 or N=3 fallback: copy nearest interior
        out[0] = out[1];
        out[N - 1] = out[N - 2];
        out[0].set_pq(out_p, out_q);
        out[N - 1].set_pq(out_p, out_q);
    }
}
/**
 * @name edthH_dmat_RSlice_inplace
 * @param in_RSlice
 * @param out_RSlice
 * @brief Compute edth_H on a single RSlice (via pointer) using the Legendre D-matrix for dz, \n
 * with safe handling of endpoints via the core method.
 */
template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthH_dmat_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t N = in_RSlice.size();
    if (N == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> df_dz_buf(N, GHPScalar<Complex>(Complex(0.0, 0.0), p, q));
    SpectralGHPVectorized::RSlice df_dz_rs(df_dz_buf.data(), N, in_RSlice.modes_, in_RSlice.r_index());

    diff_.dz_Dmatrix_RSlice(in_RSlice, df_dz_rs);
    this->edth_core_RSlice_with_extrapolation(in_RSlice, df_dz_rs, out_RSlice, EthKind::Eth);
}
/**
 * @name edthH_inplace_RSliceV
 * @param in_RSlice
 * @param out_RSlice
 * @brief Compute edth_H on a single RSlice (via spans) using the Legendre D-matrix for dz, \n
 * with safe handling of endpoints via the core method.
 */
template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthH_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> df_dz(Nz, GHPScalar<Complex>(Complex(0.0, 0.0), p, q));
    SpectralGHPVectorized::RSlice df_dz_slice(df_dz.data(), Nz, in_RSlice.modes_, in_RSlice.r_index());

    std::span<const GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> df_dz_span(df_dz_slice.data_ptr, Nz);

    diff_.dz_Dmatrix(in_span, df_dz_span);

    this->edth_core_RSlice_with_extrapolation(in_RSlice, df_dz_slice, out_RSlice, EthKind::Eth);
}

// -----------------------------------------------------------------------------
// edthBarH using dz_Dmatrix on spans (pole-safe via core)
// -----------------------------------------------------------------------------
template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthBarH_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> df_dz(Nz, GHPScalar<Complex>(Complex(0.0, 0.0), p, q));
    SpectralGHPVectorized::RSlice df_dz_slice(df_dz.data(), Nz, in_RSlice.modes_, in_RSlice.r_index());

    std::span<const GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> df_dz_span(df_dz_slice.data_ptr, Nz);

    diff_.dz_Dmatrix(in_span, df_dz_span);

    this->edth_core_RSlice_with_extrapolation(in_RSlice, df_dz_slice, out_RSlice, EthKind::EthBar);
}

/**
 * @name edthH_bary_inplace_RSliceV
 * @param in_RSlice
 * @param out_RSlice
 * @brief Compute edth_H on a single RSlice using barycentric dz, \n
 * with safe handling of endpoints via the core method. \n
 */
template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthH_bary_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& f,
        SpectralGHPVectorized::RSlice& out) const
{
    assert(f.size() == out.size());
    const size_t N = f.size();
    if (N == 0) return;

    const int p = f[0].p();
    const int q = f[0].q();

    std::vector<GHPScalar<Complex>> df_dz_buf(N, GHPScalar<Complex>(Complex(0.0, 0.0), p, q));
    SpectralGHPVectorized::RSlice df_dz_rs(df_dz_buf.data(), N, f.modes_, f.r_index());

    // Make sure diff_.z_weights() are stable weights for your LGL nodes.
    diff_.dz_barycentric_RSlice(f, df_dz_rs);

    this->edth_core_RSlice_with_extrapolation(f, df_dz_rs, out, EthKind::Eth);
}


/**
 * @name thorn_inplace_ZSliceV
 *
 **/
template <>
void KinnersleyHeldOperators<OutgoingCoords>::thorn_inplace_ZSliceV(
        const SpectralGHPVectorized::ZSlice& in_ZSlice,
        SpectralGHPVectorized::ZSlice& out_ZSlice) const
{
    assert(in_ZSlice.size() == out_ZSlice.size());
    const size_t Nr = in_ZSlice.size();
    if (Nr == 0) return;

    const int p = in_ZSlice[0].p();
    const int q = in_ZSlice[0].q();
    const size_t iz = in_ZSlice.z_index();

    std::vector<GHPScalar<Complex>> dr_buf(
            Nr, GHPScalar<Complex>(teuk::zeroC, p, q)
    );

    SpectralGHPVectorized::ZSlice dr_ZSlice(
            dr_buf.data(),
            Nr,
            1,
            in_ZSlice.modes_,
            iz,
            in_ZSlice.has_omega_mk() ? in_ZSlice.omega_mk() : Real(0),
            in_ZSlice.has_omega_mk()
    );

    diff_.dr_Dmatrix_ZSlice(in_ZSlice, dr_ZSlice);

    for (size_t ir = 0; ir < Nr; ++ir) {
        // Example placeholder:
        const auto val = dr_ZSlice[ir].value();

        out_ZSlice[ir] = GHPScalar<Complex>(val, p+1, q+1);
    }
}


/**
 * @name thornPHr_inplace_RSliceV
 * @param in_RSlice
 * @param out_RSlice
 *
 **/
template <>
void KinnersleyHeldOperators<OutgoingCoords>::thornPHr_inplace_RSliceV(
        const KinnersleyTetrad<OutgoingCoords>& ktet,
        const SpectralGHPVectorized::RSlice& in_RSlice,
        const SpectralGHPVectorized::RSlice& dr_in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    assert(dr_in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const auto& met = ktet.get_metric();
    // create lambda function for Delta(r,z)
    const size_t r_idx = in_RSlice.r_index();
    const Real rval = diff_.cl_nodes()[r_idx]; // Chebyshev nodes in r, but we want the actual r value at this index
    const Real Delta = met.Delta(rval);

    // pull omega data from slice metadata and convert to -i*omega factor for time derivative in Fourier domain
    const Complex iomega = Complex(0, 1) * in_RSlice.omega_mk();

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

#pragma omp parallel for default(none) shared(in_RSlice, met, rval, Delta, dr_in_RSlice, out_RSlice) firstprivate(Nz, iomega, p, q)
    for (size_t i = 0; i < Nz; ++i) {
        const Real zval = diff_.lgl_nodes()[i];
        const Real Sigma = met.Sigma(rval, zval);
        const Real preDr = -Delta/(2.0_r*Sigma);
        out_RSlice[i].value() = -iomega * in_RSlice[i].value() + preDr*dr_in_RSlice[i].value();
        out_RSlice[i].set_pq(p - 1, q - 1);
    }
}

template <>
void KinnersleyHeldOperators<OutgoingCoords>::thornPH_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice& in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    // pull omega data from slice metadata and convert to -i*omega factor for time derivative in Fourier domain
    const Complex iomega = Complex(0, 1) * in_RSlice.omega_mk();

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();

#pragma omp parallel for default(none) shared(in_RSlice, out_RSlice) firstprivate(Nz, iomega, p, q)
    for (size_t i = 0; i < Nz; ++i) {
        out_RSlice[i].value() = -iomega * in_RSlice[i].value();
        out_RSlice[i].set_pq(p - 1, q - 1);
    }
}