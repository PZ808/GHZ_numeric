//
// Created by Peter Zimmerman on 04.03.26.
//

#include "../include/KinnersleySpectralHeldOperators.hpp"
#include "../include/GHPScalars.hpp"
#include "../sand/SpectralGHPField.hpp"
#include "../include/SpectralDiffer.hpp"

// -----------------------------------------------------------------------------
// In-place edthH using Legendre D-matrix derivative
// -----------------------------------------------------------------------------
template<> void KinnersleyHeldOperators<OutgoingCoords>::edthH_dmat_RSlice_inplace(
        const SpectralGHPVectorized::RSlice &in_RSlice,
        SpectralGHPVectorized::RSlice &out_RSlice) const
{
    // safety
    assert(in_RSlice.size() == in_RSlice.size());
    const size_t N = in_RSlice.size();
    if (N == 0) return;

    // metadata
    const Real a = kin_tetrad_.a();
    const int p = in_RSlice[0].p();
    const int q = in_RSlice[0].q();
    const Complex s = Complex((p - q) * 0.5, 0.0); // spin weight
    const Real m = Real(in_RSlice.m());
    const Real omega = in_RSlice.omega_mk();
    const Real aw = a*omega;

    // Temporary buffer for df/dz
    std::vector<GHPScalar<Complex>> df_dz_buf;
    df_dz_buf.assign(N, GHPScalar<Complex>(Complex(0.0, 0.0), p, q));

    SpectralGHPVectorized::RSlice df_dz_rs(df_dz_buf.data(),
                                           N, in_RSlice.modes_, in_RSlice.r());

    // compute derivative using D-matrix into df_dz_buf
    diff_.dz_Dmatrix_RSlice(in_RSlice, df_dz_rs); // assumes this signature exists

    // compute edth_H into out (spin weights p+0, q-2)
    const Complex pref = Complex(-1.0 / std::sqrt(2.0), 0.0);

    for (int i = 0; i < N; ++i) {
        const Real z = diff_.lgl_nodes()[i];
        const Real factor = std::sqrt(std::max<Real>(0.0, 1.0 - z*z));

        const Complex fval = in_RSlice[i].value();
        const Complex dfz  = df_dz_rs[i].value();
        const Complex dfphi = Complex(0.0, 1.0) * Real(m) * fval;

        Complex val = pref * (
                -factor * dfz
                + Complex(0.0, 1.0) * dfphi / factor
                - s * z * fval / factor
                + Complex(aw * factor, 0.0) * fval
        );

        out_RSlice[i] = GHPScalar<Complex>(val, p+0, q - 2);
    }
}

template<> void KinnersleyHeldOperators<OutgoingCoords>::edthH_inplace_RSliceV(
        const spectral::SpectralGHPVectorized::RSlice &in_RSlice,
        spectral::SpectralGHPVectorized::RSlice &out_RSlice) const
{
    // sanity check: sizes must match
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();

    // extract GHP weights and mode info
    int p = in_RSlice[0].p();
    int q = in_RSlice[0].q();
    Complex s = Complex((p - q) / 2., 0);
    Real m = Real(in_RSlice.m());
    const Real omega = in_RSlice.omega_mk();
    const Real aw = a_*omega;
    // temporary derivative along z
    std::vector<GHPScalar<Complex>> df_dz(Nz);
    // wrap as RSlice at the same r
    spectral::SpectralGHPVectorized::RSlice df_dz_slice(df_dz.data(), Nz,
                                                        in_RSlice.modes_,
                                                        in_RSlice.r());

    // compute dz derivative
    std::span<GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> df_dz_span(df_dz_slice.data_ptr, Nz);

    diff_.dz_Dmatrix(in_span, df_dz_span);  // computes derivative along z into df_dz_span

    // main edthH computation
#pragma omp parallel for default(none) shared(Nz,p, q, out_RSlice, in_RSlice, df_dz_slice, m, s, aw)
    for (size_t i = 0; i < Nz; ++i)
    {
        Real z = diff_.lgl_nodes()[i];                     // Legendre-Gauss-Lobatto nodes
        Complex factor = Complex(std::sqrt(1.0 - z * z), 0);
        Complex df_dphi = in_RSlice[i].value() * Complex(0, 1) * Real(m);

        out_RSlice[i].value() = -Complex(1.0 / std::sqrt(2.0), 0) * (
                -factor * df_dz_slice[i].value()
                + df_dphi / factor * Complex(0, 1)
                - s * z * in_RSlice[i].value() / factor
                + aw * factor * in_RSlice[i].value()
        );

        // update GHP weights
        out_RSlice[i].set_pq(p + 0, q - 2);
    }
}

template<> void KinnersleyHeldOperators<OutgoingCoords>::edthBarH_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice &in_RSlice,
        SpectralGHPVectorized::RSlice &out_RSlice) const
{
    // sanity check: sizes must match
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();

    // extract GHP weights and mode info
    int p = in_RSlice[0].p();
    int q = in_RSlice[0].q();
    Complex s = Complex((p - q) / 2., 0);
    Real m = Real(in_RSlice.m());
    const Real omega = in_RSlice.omega_mk();
    const Real aw = a_*omega;

    // temporary derivative along z
    std::vector<GHPScalar<Complex>> df_dz(Nz);
    spectral::SpectralGHPVectorized::RSlice df_dz_slice(df_dz.data(), Nz,
                                                        in_RSlice.modes_,
                                                        in_RSlice.r());

    // compute dz derivative
    std::span<GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> df_dz_span(df_dz_slice.data_ptr, Nz);
    diff_.dz_Dmatrix(in_span, df_dz_span);

    // main edthBar computation
#pragma omp parallel for default(none) shared(Nz,p,q,out_RSlice,in_RSlice,df_dz_slice,m,s,aw)
    for (size_t i = 0; i < Nz; ++i)
    {
        Real z = diff_.lgl_nodes()[i];                     // Legendre-Gauss-Lobatto nodes
        Complex factor = Complex(std::sqrt(1.0 - z * z), 0);
        Complex df_dphi = in_RSlice[i].value() * Complex(0, 1) * Real(m);

        // edthBar flips the dz sign and some spin terms compared to edthH
        out_RSlice[i].value() = -Complex(1.0 / std::sqrt(2.0), 0) * (
                - factor * df_dz_slice[i].value()
                - df_dphi / factor * Complex(0, 1)
                + s * z * in_RSlice[i].value() / factor
                - aw * factor * in_RSlice[i].value()
        );

        // update GHP weights
        out_RSlice[i].set_pq(p - 2, q + 0);
    }
}


/**
 * @name thornPH_inplace_RSliceV
 *
 * @details Compute thornPH on a mode-decomposed Held quantity (independent of r) \n
 * on an RSlice using that the BL frequencies are conjugate \n
 * to both BL time and outgoing Kerr-Newman time and the simple FT identity: \n
 * thornPH f^circ = pd_usum_{nk} (f_{mk}^*e^{-iomega_{mk}u}) = -isum_{mnk}*omega_{mk} f_{mk}e^{-i*omega_{mk}u} \n
 * k = (kr,kz).
 *
 **/
template <> void KinnersleyHeldOperators<OutgoingCoords>::thornPH_inplace_RSliceV(
        const SpectralGHPVectorized::RSlice &in_RSlice,
        SpectralGHPVectorized::RSlice &out_RSlice) const
{
    // sanity check: sizes must match
    assert(in_RSlice.size() == out_RSlice.size());
    const size_t Nz = in_RSlice.size();
    const Complex iomega = Complex(0,1)*in_RSlice.omega_mk();

    // extract GHP weights and mode info
    int p = in_RSlice[0].p();
    int q = in_RSlice[0].q();

    std::vector<GHPScalar<Complex>> thornPH_mode(Nz);
    spectral::SpectralGHPVectorized::RSlice thornPH_mode_RSlice(thornPH_mode.data(), Nz,
                                                                in_RSlice.modes_,
                                                                in_RSlice.r());

    std::span<GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
    std::span<GHPScalar<Complex>> out_span(thornPH_mode_RSlice.data_ptr, Nz);

    // main thornPH computation
#pragma omp parallel for default(none) shared(Nz,p,q,out_RSlice,in_RSlice,thornPH_mode_RSlice,iomega)
    for (size_t i = 0; i < Nz; ++i)
    {
        out_RSlice[i].value() = -iomega*in_RSlice[i].value();

        out_RSlice[i].set_pq(p-1, q-1);
    }
}

// -----------------------------------------------------------------------------
// In-place edthH using barycentric dz
// -----------------------------------------------------------------------------
template <> void KinnersleyHeldOperators<OutgoingCoords>::edthH_bary_inplace(
        const SpectralGHPVectorized::RSlice &f,
        SpectralGHPVectorized::RSlice &out) const
{
    // safety
    assert(f.size() == out.size());
    const size_t N = f.size();
    if (N==0) return;

    // metadata
    const int p = f[0].p();
    const int q = f[0].q();
    const Complex s = Complex((p - q) * 0.5, 0.0);
    const Real m = Real(f.m());
    const Real aw = a_*f.omega_mk();

    // Temporary buffer for df/dz (owned here -> safe if out aliases f)
    std::vector<GHPScalar<Complex>> df_dz_buf;
    df_dz_buf.assign(N, GHPScalar<Complex>(Complex(0.0, 0.0), p, q));

    SpectralGHPVectorized::RSlice df_dz_rs(df_dz_buf.data(), N, f.modes_, f.r());

    // compute derivative in-place into df_dz_buf
    diff_.dz_barycentric_RSlice(f, df_dz_rs, diff_.z_weights()); // assumes this signature exists

    // compute edth_H into out (spin weights p+1, q-1)
    const Complex pref = Complex(-1.0 / std::sqrt(2.0), 0.0);

    for (int i = 0; i < N; ++i) {
        const Real z = diff_.lgl_nodes()[i];
        const Real factor = std::sqrt(std::max<Real>(0.0, 1.0 - z*z));

        const Complex fval = f[i].value();
        const Complex dfz  = df_dz_rs[i].value();
        const Complex dfphi = Complex(0.0, 1.0) * Real(m) * fval;

        Complex val = pref * (
                -factor * dfz
                + Complex(0.0, 1.0) * dfphi / factor
                - s * z * fval / factor
                + Complex(aw * factor, 0.0) * fval
        );

        out[i] = GHPScalar<Complex>(val, p + 0, q - 2);
    }
}


