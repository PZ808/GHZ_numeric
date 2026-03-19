//
// Created by Peter Zimmerman on 04.03.26.
//
#ifndef GHZ_SPECTRAL_KINNERSLEYHELDOPERATORS_HPP
#define GHZ_SPECTRAL_KINNERSLEYHELDOPERATORS_HPP

#pragma once
#include "ghz/core/GhzTypes.hpp"
#include "ghz/geom/KinnersleyTetrad.hpp"
#include "ghz/ghp/GHPScalars.hpp"
#include "ghz/ghp/HeldScalars.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/spectral/SpectralCoordinateMaps.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"

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
    const ghz::numeric::AffineMap1D r_map_;
    const Real a_ = kin_tetrad_.a();
public:

    KinnersleyHeldOperators(const DifferType& diff,
                            const ghp::HeldBackgroundFieldsVectorized<CoordType>& bg_helds,
                            const KinnersleyTetrad<CoordType>& kin_tetrad,
                            const ghz::numeric::AffineMap1D& r_map)
            : diff_(diff), bg_helds_(bg_helds), kin_tetrad_(kin_tetrad),
            r_map_(r_map) {};

    [[nodiscard]] const ghz::numeric::AffineMap1D& r_map() const noexcept { return r_map_; }
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

    void thornPH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                                 SpectralGHPVectorized::RSlice& out) const;

    void edthH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                               SpectralGHPVectorized::RSlice& out) const;
    void edthBarH_inplace_RSliceV(const SpectralGHPVectorized::RSlice& in,
                                  SpectralGHPVectorized::RSlice& out) const;
    void edthH_dmat_inplace_RSliceV(const SpectralGHPVectorized::RSlice &in_RSlice,
                                    SpectralGHPVectorized::RSlice &out_RSlice) const;

    void edthH_bary_inplace_RSliceV(const SpectralGHPVectorized::RSlice &f,
                                    SpectralGHPVectorized::RSlice &out) const;

    void thornPHr_inplace_RSliceV(const KinnersleyTetrad<OutgoingCoords> &ktet,
                                  const SpectralGHPVectorized::RSlice &in_RSlice,
                                  const SpectralGHPVectorized::RSlice &dr_in_RSlice,
                                  SpectralGHPVectorized::RSlice &out_RSlice) const;


    void thorn_inplace_ZSliceV(
            const SpectralGHPVectorized::ZSlice& in_ZSlice,
            SpectralGHPVectorized::ZSlice& out_ZSlice) const;

    // ConstRSlice views
    void edthH_inplace_RSliceV(const SpectralGHPVectorized::ConstRSlice& in,
                               SpectralGHPVectorized::RSlice& out) const;
    void edthBarH_inplace_RSliceV(const SpectralGHPVectorized::ConstRSlice& in,
                                  SpectralGHPVectorized::RSlice& out) const;
    void thornPH_inplace_RSliceV(const SpectralGHPVectorized::ConstRSlice& in,
                                 SpectralGHPVectorized::RSlice& out) const;
    void thornPHr_inplace_RSliceV(const KinnersleyTetrad<OutgoingCoords> &ktet,
                                  const SpectralGHPVectorized::ConstRSlice &in_RSlice,
                                  const SpectralGHPVectorized::ConstRSlice &dr_in_RSlice,
                                  SpectralGHPVectorized::RSlice &out_RSlice) const;
    // ConstZlice views
    void thornPHr_inplace_ZSliceV(const KinnersleyTetrad<OutgoingCoords> &ktet,
                                  const SpectralGHPVectorized::ConstZSlice &in_ZSlice,
                                  SpectralGHPVectorized::ZSlice &out_ZSlice) const;
    void thorn_inplace_ZSliceV(
            const SpectralGHPVectorized::ConstZSlice& in_ZSlice,
            SpectralGHPVectorized::ZSlice& out_ZSlice) const;

    // wrappers to generate complete 2d fields
    void edthH_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;
    void edthBarH_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;
    void thornPH_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;
    void thornPHr_inplace(const KinnersleyTetrad<OutgoingCoords> &ktet,
                          const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;
    void thorn_inplace(const SpectralGHPVectorized& in, SpectralGHPVectorized& out) const;

};

#endif //GHZ_SPECTRAL_KINNERSLEYHELDOPERATORS_HPP
