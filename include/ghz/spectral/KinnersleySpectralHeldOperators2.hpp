//
// Created by Peter Zimmerman on 04.03.26.
//
#ifndef HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP
#define HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP

#pragma once
#include "ghz/core/GhzTypes.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"
#include "ghz/geom/KinnersleyTetrad.hpp"

using namespace  spectral;
using namespace ghp;

namespace detail {

    static inline Complex quad_extrapolate_endpoint(
            Real z0,
            Real z1, Complex f1,
            Real z2, Complex f2,
            Real z3, Complex f3)
    {
        // Quadratic Lagrange extrapolation to z0 from points (z1,z2,z3).
        const Real d12 = (z1 - z2);
        const Real d13 = (z1 - z3);
        const Real d21 = (z2 - z1);
        const Real d23 = (z2 - z3);
        const Real d31 = (z3 - z1);
        const Real d32 = (z3 - z2);

        Complex L1 = f1 * ((z0 - z2) * (z0 - z3)) / (d12 * d13);
        Complex L2 = f2 * ((z0 - z1) * (z0 - z3)) / (d21 * d23);
        Complex L3 = f3 * ((z0 - z1) * (z0 - z2)) / (d31 * d32);
        return L1 + L2 + L3;
    }

} // namespace detail

enum class EthKind { Eth, EthBar };

template <typename CoordType, typename DifferType=spectral::SpectralDiffer>
class KinnersleyHeldOperators {

private:
    const DifferType& diff_;
    const HeldBackgroundFieldsVectorized<CoordType>& bg_helds_;
    const KinnersleyTetrad<CoordType>& kin_tetrad_;
    const Real a_ = kin_tetrad_.a();
public:

    KinnersleyHeldOperators(const DifferType& diff,
                            const ghp::HeldBackgroundFieldsVectorized<CoordType>& bg_helds,
                            const KinnersleyTetrad<CoordType>& kin_tetrad)
            : diff_(diff), bg_helds_(bg_helds), kin_tetrad_(kin_tetrad) {}

/**
 * @name eth_core_RSlice_with_extrapolation
 * @tparam CoordType coordinate type (e.g., Outgoing Boyer-Lindquist, Kerr-Schild)
 * @param in input spectral field on a single RSlice (fixed r, varying z) \n
 * @param df_dz precomputed spectral derivative d/dz of the input field on the same RSlice \n
 * @param out output spectral field where the result will be written, same size as in and df_dz \n
 * @param kind specifies whether to compute the edth or edthBar operator \n
 * @brief Compute the core part of the edth/edthBar operator on a single RSlice, with safe handling of endpoints.
 */
    template <typename InSlice>
    void edth_core_RSlice_with_extrapolation(
            const InSlice& in_RSlice,
            const spectral::SpectralGHPVectorized::RSlice& df_dz,
            spectral::SpectralGHPVectorized::RSlice& out,
            EthKind kind) const
    {
        assert(in_RSlice.size() == out.size());
        assert(in_RSlice.size() == df_dz.size());

        const size_t N=in_RSlice.size();
        if (N==0) return;

        const int p = in_RSlice[0].p();
        const int q = in_RSlice[0].q();
        const Complex s = Complex(Real(p-q)*Real(0.5),0.0);
        const Real m = Real(in_RSlice.m());
        const Real aw = a_*in_RSlice.omega_mk();

        const int out_p = (kind==EthKind::Eth) ? p : p-2;
        const int out_q = (kind==EthKind::Eth) ? q-2 : q;

        const Complex pref=Complex(-Real(1.0)/Real(std::sqrt(2.0)),0.0);

        const auto& z_nodes=diff_.lgl_nodes();
        assert(z_nodes.size()==N);

#pragma omp parallel for default(none) shared(in_RSlice,df_dz,out,z_nodes,kind) \
    firstprivate(N,s,m,aw,pref,out_p,out_q)
        for (size_t i=1; i<N-1; ++i) {
            const Real z=z_nodes[i];
            const Real fac_r=std::sqrt(std::max<Real>(Real(0.0), Real(1.0)-z*z));

            if (fac_r<=Real(0.0)) {
                out[i]=GHPScalar<Complex>(Complex(0.0,0.0),out_p,out_q);
                continue;
            }

            const Complex factor(fac_r,0.0);

            const Complex fval = in_RSlice[i].value();
            const Complex dfz = df_dz[i].value();

            Complex singular_num;
            if (kind == EthKind::Eth) {
                singular_num = Complex(-m,0.0)*fval-s*Complex(z,0.0)*fval;
            } else {
                singular_num = Complex(+m,0.0)*fval+s*Complex(z,0.0)*fval;
            }
            const Complex singular = singular_num/factor;

            const Complex aw_term = (kind==EthKind::Eth)
                                  ? Complex(aw,0.0)*factor*fval
                                  : Complex(-aw,0.0)*factor*fval;

            const Complex dz_term = -factor*dfz;

            out[i]=GHPScalar<Complex>(pref*(dz_term+singular+aw_term),out_p,out_q);
        }

        if (N==1) {
            out[0]=GHPScalar<Complex>(Complex(0.0,0.0),out_p,out_q);
            return;
        }

        if (N>=4) {
            out[0].value()=detail::quad_extrapolate_endpoint(
                    z_nodes[0],
                    z_nodes[1],out[1].value(),
                    z_nodes[2],out[2].value(),
                    z_nodes[3],out[3].value());
            out[0].set_pq(out_p,out_q);

            out[N-1].value()=detail::quad_extrapolate_endpoint(
                    z_nodes[N-1],
                    z_nodes[N-2],out[N-2].value(),
                    z_nodes[N-3],out[N-3].value(),
                    z_nodes[N-4],out[N-4].value());
            out[N-1].set_pq(out_p,out_q);
        } else {
            out[0]=out[1];
            out[N-1]=out[N-2];
            out[0].set_pq(out_p,out_q);
            out[N-1].set_pq(out_p,out_q);
        }
    }
    //void edth_core_RSlice_with_extrapolation(const spectral::SpectralGHPVectorized::RSlice& in, const spectral::SpectralGHPVectorized::RSlice& df_dz, spectral::SpectralGHPVectorized::RSlice& out, EthKind kind) const;
    template <typename InSlice>
    void edthH_inplace_RSliceV(const InSlice& in_RSlice,
                               SpectralGHPVectorized::RSlice& out_RSlice) const;

    template <typename InSlice>
    void edthBarH_inplace_RSliceV(const InSlice& in_RSlice,
                                  SpectralGHPVectorized::RSlice& out_RSlice) const;

    template <typename InSlice>
    void thornPH_inplace_RSliceV(const InSlice& in_RSlice,
                                 SpectralGHPVectorized::RSlice& out_RSlice) const;

    template <typename InSlice>
    void thorn_inplace_ZSliceV(const InSlice& in_ZSlice,
                               SpectralGHPVectorized::ZSlice& out_ZSlice) const;

    template <typename InSlice, typename DrSlice>
    void thornPHr_inplace_RSliceV(const KinnersleyTetrad<OutgoingCoords>& ktet,
                                  const InSlice& in_RSlice,
                                  const DrSlice& dr_in_RSlice,
                                  SpectralGHPVectorized::RSlice& out_RSlice) const;

    template <typename InSlice>
    void thornPHr_inplace_ZSliceV(const KinnersleyTetrad<OutgoingCoords>& ktet,
                                  const InSlice& in_ZSlice,
                                  SpectralGHPVectorized::ZSlice& out_ZSlice) const;

    void edthH_dmat_inplace_RSliceV(const SpectralGHPVectorized::RSlice &in_RSlice,
                                    SpectralGHPVectorized::RSlice &out_RSlice) const;

    void edthH_bary_inplace_RSliceV(const SpectralGHPVectorized::RSlice &f,
                                    SpectralGHPVectorized::RSlice &out) const;
    template <typename InSlice>
    void edthH_inplace_RSliceV(const InSlice& in_RSlice,
                               SpectralGHPVectorized::RSlice& out_RSlice) const
    {
        assert(in_RSlice.size()==out_RSlice.size());
        const size_t Nz = in_RSlice.size();
        if (Nz==0) return;

        const int p = in_RSlice[0].p();
        const int q = in_RSlice[0].q();

        std::vector<GHPScalar<Complex>> df_dz(
                Nz, GHPScalar<Complex>(Complex(0.0, 0.0), p, q)
        );

        SpectralGHPVectorized::RSlice df_dz_slice(
                df_dz.data(),
                Nz,
                in_RSlice.modes_,
                in_RSlice.r_index(),
                in_RSlice.has_omega_mk() ? in_RSlice.omega_mk() : Real(0),
                in_RSlice.has_omega_mk()
        );

        std::span<const GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
        std::span<GHPScalar<Complex>> df_dz_span(df_dz_slice.data_ptr, Nz);

        diff_.dz_Dmatrix(in_span, df_dz_span);

        this->edth_core_RSlice_with_extrapolation(in_RSlice, df_dz_slice, out_RSlice, EthKind::Eth);
    }
    template <typename InSlice>
    void edthBarH_inplace_RSliceV(const InSlice& in_RSlice,
                                  SpectralGHPVectorized::RSlice& out_RSlice) const
    {
        assert(in_RSlice.size()==out_RSlice.size());
        const size_t Nz = in_RSlice.size();
        if (Nz==0) return;

        const int p = in_RSlice[0].p();
        const int q = in_RSlice[0].q();

        std::vector<GHPScalar<Complex>> df_dz(
                Nz, GHPScalar<Complex>(Complex(0.0, 0.0), p, q)
        );

        SpectralGHPVectorized::RSlice df_dz_slice(
                df_dz.data(),
                Nz,
                in_RSlice.modes_,
                in_RSlice.r_index(),
                in_RSlice.has_omega_mk() ? in_RSlice.omega_mk() : Real(0),
                in_RSlice.has_omega_mk()
        );

        std::span<const GHPScalar<Complex>> in_span(in_RSlice.data_ptr, Nz);
        std::span<GHPScalar<Complex>> df_dz_span(df_dz_slice.data_ptr, Nz);

        diff_.dz_Dmatrix(in_span, df_dz_span);

        this->edth_core_RSlice_with_extrapolation(in_RSlice, df_dz_slice, out_RSlice, EthKind::EthBar);
    }
    template <typename InSlice>
    void thornPH_inplace_RSliceV(const InSlice& in_RSlice,
                                 SpectralGHPVectorized::RSlice& out_RSlice) const
    {
        assert(in_RSlice.size()==out_RSlice.size());
        const size_t Nz = in_RSlice.size();
        if (Nz==0) return;

        const Complex iomega = Complex(0, 1)*in_RSlice.omega_mk();

        const int p = in_RSlice[0].p();
        const int q = in_RSlice[0].q();

#pragma omp parallel for default(none) shared(in_RSlice, out_RSlice) firstprivate(Nz, iomega, p, q)
        for (size_t i = 0; i<Nz; ++i) {
            out_RSlice[i].value() = -iomega*in_RSlice[i].value();
            out_RSlice[i].set_pq(p-1, q-1);
        }
    }
    template <typename InSlice>
    void thorn_inplace_ZSliceV(const InSlice& in_ZSlice,
                               SpectralGHPVectorized::ZSlice& out_ZSlice) const
    {
        assert(in_ZSlice.size()==out_ZSlice.size());
        const size_t Nr = in_ZSlice.size();
        if (Nr==0) return;

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

        for (size_t ir = 0; ir<Nr; ++ir) {
            out_ZSlice[ir] = GHPScalar<Complex>(dr_ZSlice[ir].value(), p+1, q+1);
        }
    }

    // wrappers to generate complete 2d fields
    void edthH_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;
    void edthBarH_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;
    void ThornPH_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;
    void ThornPHr_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;
    void Thorn_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;

};

#endif //HS1DATA_DATA_KINNERSLEYHELDOPERATORS_HPP
