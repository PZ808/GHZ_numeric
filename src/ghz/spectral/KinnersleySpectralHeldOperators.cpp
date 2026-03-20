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
// ConstRSlice version of edthBarH_inplace_RSliceV using dz_Dmatrix on spans (pole-safe via core)
template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthH_inplace_RSliceV(
        const SpectralGHPVectorized::ConstRSlice& in_RSlice,
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
// -----------------------------------------------------------------------------
// edthBarH using dz_Dmatrix on spans (pole-safe via core)
// -----------------------------------------------------------------------------
template <>
void KinnersleyHeldOperators<OutgoingCoords>::edthBarH_inplace_RSliceV(
        const SpectralGHPVectorized::ConstRSlice& in_RSlice,
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

    std::vector<GHPScalar<Complex>> dr_buf( Nr, GHPScalar<Complex>(teuk::zeroC, p, q) );

    SpectralGHPVectorized::ZSlice dr_ZSlice(
            dr_buf.data(),
            Nr,
            1,
            in_ZSlice.modes_,
            iz,
            in_ZSlice.has_omega_mk() ? in_ZSlice.omega_mk() : Real(0),
            in_ZSlice.has_omega_mk() );

    diff_.dx_Dmatrix_ZSlice(in_ZSlice, dr_ZSlice);
    const Real dxdr = r_map_.dxdr();

    for (size_t ir = 0; ir < Nr; ++ir) {
        const auto val = dxdr*dr_ZSlice[ir].value();
        out_ZSlice[ir] = GHPScalar<Complex>(val, p+1, q+1);
    }
}
// const version
template <>
void KinnersleyHeldOperators<OutgoingCoords>::thorn_inplace_ZSliceV(
        const SpectralGHPVectorized::ConstZSlice& in_ZSlice,
        SpectralGHPVectorized::ZSlice& out_ZSlice) const
{
    assert(in_ZSlice.size() == out_ZSlice.size());
    const size_t Nr = in_ZSlice.size();
    if (Nr == 0) return;

    const int p = in_ZSlice[0].p();
    const int q = in_ZSlice[0].q();
    const size_t iz = in_ZSlice.z_index();

    std::vector<GHPScalar<Complex>> dr_buf( Nr, GHPScalar<Complex>(teuk::zeroC, p, q) );

    SpectralGHPVectorized::ZSlice dr_ZSlice(
            dr_buf.data(),
            Nr,
            1,
            in_ZSlice.modes_,
            iz,
            in_ZSlice.has_omega_mk() ? in_ZSlice.omega_mk() : Real(0),
            in_ZSlice.has_omega_mk() );

    diff_.dx_Dmatrix_ZSlice(in_ZSlice, dr_ZSlice);
    const Real dxdr = r_map_.dxdr();

    for (size_t ir = 0; ir < Nr; ++ir) {
        // apply jacobian to change from x to r
        const auto val = dxdr * dr_ZSlice[ir].value();

        out_ZSlice[ir] = GHPScalar<Complex>(val, p+1, q+1);
    }
}


/**
 * @name thornPHr_inplace_RSliceV
 * @param in_RSlice
 * @param out_RSlice
 * @brief Compute thornPHr on a single RSlice (fixed r val) using dr_Dmatrix for the radial derivative, \n
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
    const Real x = diff_.cl_nodes()[r_idx]; // Chebyshev nodes in x, but we want the actual r value at this index
    const Real rval = r_map_.toPhysical(x);
    const Real Delta = met.Delta(rval);

    // pull omega data from slice metadata and convert to -i*omega factor for time derivative in Fourier domain
    const Complex iomega = Complex(0, 1) * in_RSlice.omega_mk();

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();
    const Real dxdr = r_map_.dxdr();

#pragma omp parallel for default(none) shared(in_RSlice, met, dxdr, rval, Delta, dr_in_RSlice, out_RSlice) firstprivate(Nz, iomega, p, q)
    for (size_t i = 0; i < Nz; ++i) {
        const Real zval = diff_.lgl_nodes()[i];
        const Real Sigma = met.Sigma(rval, zval);
        const Real preDr = -Delta/(2.0_r*Sigma);
        out_RSlice[i].value() = -iomega * in_RSlice[i].value() + preDr*dxdr*dr_in_RSlice[i].value();
        out_RSlice[i].set_pq(p - 1, q - 1);
    }
}

template <>
void KinnersleyHeldOperators<OutgoingCoords>::thornPHr_inplace_RSliceV(
        const KinnersleyTetrad<OutgoingCoords>& ktet,
        const SpectralGHPVectorized::ConstRSlice& in_RSlice,
        const SpectralGHPVectorized::ConstRSlice& dr_in_RSlice,
        SpectralGHPVectorized::RSlice& out_RSlice) const
{
    assert(in_RSlice.size() == out_RSlice.size());
    assert(dr_in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    if (Nz == 0) return;

    const auto& met = ktet.get_metric();
    // create lambda function for Delta(r,z)
    const size_t r_idx = in_RSlice.r_index();
    const Real x = diff_.cl_nodes()[r_idx]; // Chebyshev nodes in r, but we want the actual r value at this index
    const Real rval = r_map_.toPhysical(x);
    const Real Delta = met.Delta(rval);

    // pull omega data from slice metadata and convert to -i*omega factor for time derivative in Fourier domain
    const Complex iomega = Complex(0, 1) * in_RSlice.omega_mk();

    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();
    const Real dxdr = r_map_.dxdr();

#pragma omp parallel for default(none) shared(in_RSlice, met, rval, Delta, dr_in_RSlice, dxdr, out_RSlice) firstprivate(Nz, iomega, p, q)
    for (size_t i = 0; i < Nz; ++i) {
        const Real zval = diff_.lgl_nodes()[i];
        const Real Sigma = met.Sigma(rval, zval);
        const Real preDr = -Delta/(2.0_r*Sigma);
        out_RSlice[i].value() = -iomega * in_RSlice[i].value() + preDr*dxdr*dr_in_RSlice[i].value();
        out_RSlice[i].set_pq(p - 1, q - 1);
    }
}
template <>
void KinnersleyHeldOperators<OutgoingCoords>::thornPHr_inplace_ZSliceV(
        const KinnersleyTetrad<OutgoingCoords>& ktet,
        const SpectralGHPVectorized::ConstZSlice& in_ZSlice,
        SpectralGHPVectorized::ZSlice& out_ZSlice) const
{
    assert(in_ZSlice.size() == out_ZSlice.size());
    const size_t Nr = in_ZSlice.size();
    const auto& met = ktet.get_metric();
    // create lambda function for Delta(r,z)
    const size_t z_idx = in_ZSlice.z_index();
    const Real zval = diff_.lgl_nodes()[z_idx]; // Chebyshev nodes in r, but we want the actual r value at this index


    // pull omega data from slice metadata and convert to -i*omega factor for time derivative in Fourier domain
    const Complex iomega = Complex(0, 1) * in_ZSlice.omega_mk();

    const int p = in_ZSlice[0].p();
    const int q = in_ZSlice[0].q();
    std::vector<GHPScalar<Complex>> dr_buf( Nr, GHPScalar<Complex>(teuk::zeroC, p, q) );

    SpectralGHPVectorized::ZSlice dr_ZSlice(
            dr_buf.data(),
            Nr,
            1,
            in_ZSlice.modes_,
            z_idx,
            in_ZSlice.has_omega_mk() ? in_ZSlice.omega_mk() : Real(0),
            in_ZSlice.has_omega_mk() );

    diff_.dx_Dmatrix_ZSlice(in_ZSlice, dr_ZSlice);
    const Real dxdr = r_map_.dxdr();

//#pragma omp parallel for default(none) shared(in_ZSlice, met, zval, dr_ZSlice, out_ZSlice) firstprivate(Nr, iomega, p, q)
    for (size_t i = 0; i < Nr; ++i) {
        const Real x = diff_.cl_nodes()[i];
        const Real rval = r_map_.toPhysical(x);
        const Real Delta = met.Delta(rval);
        const Real Sigma = met.Sigma(rval, zval);
        const Real preDr = -Delta/(2.0_r*Sigma);
        const auto dfdr = dxdr*dr_ZSlice[i].value();
        out_ZSlice[i].value() = -iomega * in_ZSlice[i].value() + preDr*dfdr;
        out_ZSlice[i].set_pq(p - 1, q - 1);
    }
}

/**
 * @name thornPH_inplace_RSliceV
 * @brief   strict Held thornPH operator on R-slices, which reduces to -i*omega in Fourier domain for time derivatives, \n
 *
 */
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
template <>
void KinnersleyHeldOperators<OutgoingCoords>::thornPH_inplace_RSliceV(
        const SpectralGHPVectorized::ConstRSlice& in_RSlice,
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

//-----------------------------------------------------------------------------
//  wrappers
// -----------------------------------------------------------------------------
// wrapper to apply edthH to all R-slices of a full SpectralGHPVectorized field,
// using the RSliceV version above
template <>
void KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::edthH_inplace(
        const SpectralGHPVectorized& in,
        SpectralGHPVectorized& out) const {
    assert(in.Nr() == out.Nr() && in.Nz() == out.Nz());

    for (size_t ir = 0; ir < in.Nr(); ++ir) {
        auto in_r  = in.slice_R(ir);   // const overload -> ConstRSlice
        auto out_r = out.slice_R(ir);  // non-const overload -> RSlice
        edthH_inplace_RSliceV(in_r, out_r);
    }
}

// wrapper to apply edthBarH to all R-slices of a full SpectralGHPVectorized field,
// using the RSliceV version above
template <>
void KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::edthBarH_inplace(
        const SpectralGHPVectorized& in,
        SpectralGHPVectorized& out) const {
    assert(in.Nr() == out.Nr() && in.Nz() == out.Nz());

    for (size_t ir = 0; ir < in.Nr(); ++ir) {
        auto in_r  = in.slice_R(ir);
        auto out_r = out.slice_R(ir);
        edthBarH_inplace_RSliceV(in_r, out_r);
    }
}
template <>
void  KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::thornPH_inplace(
        const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const {
    assert(in.Nr() == out.Nr() && in.Nz() == out.Nz());

    for (size_t ir = 0; ir < in.Nr(); ++ir) {
        auto in_r  = in.slice_R(ir);
        auto out_r = out.slice_R(ir);
        thornPH_inplace_RSliceV(in_r, out_r);
    }
}
template <>
void  KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::thornPHr_inplace(
        const KinnersleyTetrad<OutgoingCoords> &ktet,
        const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const {
    assert(in.Nr() == out.Nr() && in.Nz() == out.Nz());

    for (size_t iz = 0; iz < in.Nz(); ++iz) {
        auto in_z  = in.slice_Z(iz);
        auto out_z = out.slice_Z(iz);
        thornPHr_inplace_ZSliceV(ktet, in_z, out_z);
    }
}
template <>
void  KinnersleyHeldOperators<OutgoingCoords, spectral::SpectralDiffer>::thorn_inplace(
        const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const {
    assert(in.Nr() == out.Nr() && in.Nz() == out.Nz());

    for (size_t iz = 0; iz < in.Nz(); ++iz) {
        auto in_z  = in.slice_Z(iz);
        auto out_z = out.slice_Z(iz);
        thorn_inplace_ZSliceV(in_z, out_z);
    }
}