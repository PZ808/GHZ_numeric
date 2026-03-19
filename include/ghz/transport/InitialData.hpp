//
// Created by Peter Zimmerman on 13.03.26.
//

#ifndef GHZ_NUMERIC_INITIALDATA_HPP
#define GHZ_NUMERIC_INITIALDATA_HPP
#pragma once

#include "ghz/core/GhzTypes.hpp"
#include "ghz/spectral/SpectralGHPFieldVectorized.hpp"
#include "ghz/spectral/SpectralDiffer.hpp"

#include <vector>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace ghz {


    class InitialDataXmmbar {
    public:
        using Real     = teuk::Real;
        using Complex  = teuk::Complex;
        using StateVec = std::vector<Real>;

        struct BoundaryData {
            // all on the z-grid at r = r_min
            std::span<const Real> Tll_deltaPrime_minus;   // X_mmbar^-(z)
            std::span<const Real> Tll_delta_minus;        // part of d_r X_mmbar^-
            std::span<const Real> rhoPlusRhob_at_rmin;    // (rho + rhob)|_{rmin}
        };

        explicit InitialDataXmmbar(std::vector<Real> zvals, std::size_t ir_min = 0)
                : zvals_(std::move(zvals)), ir_min_(ir_min)
        {
            if (zvals_.empty()) {
                throw std::runtime_error("InitialDataXmmbar: zvals cannot be empty.");
            }
        }

        [[nodiscard]] std::size_t size() const noexcept { return zvals_.size(); }
        [[nodiscard]] std::size_t ir_min() const noexcept { return ir_min_; }
        [[nodiscard]] const std::vector<Real>& zvals() const noexcept { return zvals_; }

        [[nodiscard]] Real X0(std::size_t iz, const BoundaryData& bd) const {
            check_sizes_(bd);
            return bd.Tll_deltaPrime_minus[iz];
        }

        [[nodiscard]] Real dX0(std::size_t iz, const BoundaryData& bd) const {
            check_sizes_(bd);
            return bd.Tll_delta_minus[iz]
                   + bd.rhoPlusRhob_at_rmin[iz] * bd.Tll_deltaPrime_minus[iz];
        }

        [[nodiscard]] StateVec make_state(std::size_t iz, const BoundaryData& bd) const {
            return StateVec{ X0(iz, bd), dX0(iz, bd) };
        }

        [[nodiscard]] std::vector<StateVec> make_all_states(const BoundaryData& bd) const {
            check_sizes_(bd);
            std::vector<StateVec> out;
            out.reserve(zvals_.size());

            for (std::size_t iz = 0; iz < zvals_.size(); ++iz) {
                out.push_back(make_state(iz, bd));
            }
            return out;
        }

        template <typename FieldX, typename FieldDX>
        void fill_boundary_row(FieldX& Xmmbar,
                               FieldDX& dXmmbar_dr,
                               const BoundaryData& bd) const
        {
            check_sizes_(bd);
            check_field_shape_(Xmmbar, "Xmmbar");
            check_field_shape_(dXmmbar_dr, "dXmmbar_dr");

            for (std::size_t iz = 0; iz < zvals_.size(); ++iz) {
                Xmmbar(ir_min_, iz)     = X0(iz, bd);
                dXmmbar_dr(ir_min_, iz) = dX0(iz, bd);
            }
        }

    private:
        std::vector<Real> zvals_;
        std::size_t ir_min_{0};

        void check_sizes_(const BoundaryData& bd) const {
            const std::size_t Nz = zvals_.size();
            if (bd.Tll_deltaPrime_minus.size() != Nz ||
                bd.Tll_delta_minus.size()      != Nz ||
                bd.rhoPlusRhob_at_rmin.size()  != Nz)
            {
                throw std::runtime_error(
                        "InitialDataXmmbar: boundary arrays must all have size zvals.size()."
                );
            }
        }

        template <typename Field>
        void check_field_shape_(const Field& field, const std::string& name) const {
            if (ir_min_ >= field.Nr()) {
                throw std::runtime_error("InitialDataXmmbar: ir_min out of range for " + name);
            }
            if (field.Nz() != zvals_.size()) {
                throw std::runtime_error("InitialDataXmmbar: Nz mismatch for " + name);
            }
        }
    };



    class InitialDataXnm {
    public:
        using Real     = teuk::Real;
        using Complex  = teuk::Complex;
        using StateVec = std::vector<Real>;

        struct BoundaryData {
            // all on the z-grid at r = r_min
            std::span<const Complex> Tlm_deltaPrime_minus;
            std::span<const Complex> Tlm_delta_minus;
            std::span<const Complex> tedth_Tlm_deltaPrime_minus;
            std::span<const Complex> rho_at_rmin;
        };

        explicit InitialDataXnm(std::vector<Real> zvals, std::size_t ir_min = 0)
                : zvals_(std::move(zvals)), ir_min_(ir_min)
        {
            if (zvals_.empty()) {
                throw std::runtime_error("InitialDataXnm: zvals cannot be empty.");
            }
        }

        [[nodiscard]] std::size_t size() const noexcept { return zvals_.size(); }
        [[nodiscard]] std::size_t ir_min() const noexcept { return ir_min_; }
        [[nodiscard]] const std::vector<Real>& zvals() const noexcept { return zvals_; }

        [[nodiscard]] Complex X0(std::size_t iz, const BoundaryData& bd) const {
            check_sizes_(bd);
            return Real(2) * bd.Tlm_deltaPrime_minus[iz];
        }

        [[nodiscard]] Complex dX0(std::size_t iz, const BoundaryData& bd) const {
            check_sizes_(bd);

            const Complex rho  = bd.rho_at_rmin[iz];
            const Complex rhob = std::conj(rho);

            return Real(4) * rho  * bd.Tlm_deltaPrime_minus[iz]
                   + Real(2) * bd.Tlm_delta_minus[iz]
                   - rhob * bd.tedth_Tlm_deltaPrime_minus[iz];
        }

        [[nodiscard]] StateVec make_state(std::size_t iz, const BoundaryData& bd) const {
            const Complex x0  = X0(iz, bd);
            const Complex dx0 = dX0(iz, bd);

            return StateVec{
                    x0.real(),
                    dx0.real(),
                    x0.imag(),
                    dx0.imag()
            };
        }

        [[nodiscard]] std::vector<StateVec> make_all_states(const BoundaryData& bd) const {
            check_sizes_(bd);
            std::vector<StateVec> out;
            out.reserve(zvals_.size());

            for (std::size_t iz = 0; iz < zvals_.size(); ++iz) {
                out.push_back(make_state(iz, bd));
            }
            return out;
        }

        template <typename FieldX, typename FieldDX>
        void fill_boundary_row(FieldX& Xnm,
                               FieldDX& dXnm_dr,
                               const BoundaryData& bd) const
        {
            check_sizes_(bd);
            check_field_shape_(Xnm, "Xnm");
            check_field_shape_(dXnm_dr, "dXnm_dr");

            for (std::size_t iz = 0; iz < zvals_.size(); ++iz) {
                Xnm(ir_min_, iz)    = X0(iz, bd);
                dXnm_dr(ir_min_, iz) = dX0(iz, bd);
            }
        }

    private:
        std::vector<Real> zvals_;
        std::size_t ir_min_{0};

        void check_sizes_(const BoundaryData& bd) const {
            const std::size_t Nz = zvals_.size();
            if (bd.Tlm_deltaPrime_minus.size()       != Nz ||
                bd.Tlm_delta_minus.size()            != Nz ||
                bd.tedth_Tlm_deltaPrime_minus.size() != Nz ||
                bd.rho_at_rmin.size()                != Nz)
            {
                throw std::runtime_error(
                        "InitialDataXnm: boundary arrays must all have size zvals.size()."
                );
            }
        }

        template <typename Field>
        void check_field_shape_(const Field& field, const std::string& name) const {
            if (ir_min_ >= field.Nr()) {
                throw std::runtime_error("InitialDataXnm: ir_min out of range for " + name);
            }
            if (field.Nz() != zvals_.size()) {
                throw std::runtime_error("InitialDataXnm: Nz mismatch for " + name);
            }
        }
    };


    class InitialDataXnn {
    public:
        using Real     = teuk::Real;
        using Complex  = teuk::Complex;
        using StateVec = std::vector<Real>;

        struct BoundaryData {
            // all quantities evaluated on the z-grid at r = r_min

            std::span<const Complex> Tln_delta_minus;                 // T_ln^{(δ,-)}
            std::span<const Complex> Tll_delta_minus;                 // T_ll^{(δ,-)}
            std::span<const Complex> du_Tll_deltaPrime_minus;         // ∂_u T_ll^{(δ',-)}
            std::span<const Complex> Tll_deltaPrime_minus;            // T_ll^{(δ',-)}

            std::span<const Complex> tedth_Tlm_deltaPrime_minus;      // \tilde{edth} T_lm^{(δ',-)}
            std::span<const Complex> tedth_Tlmbar_deltaPrime_minus;   // \tilde{edth} T_l\bar{m}^{(δ',-)}

            std::span<const Complex> Tlm_deltaPrime_minus;            // T_lm^{(δ',-)}
            std::span<const Complex> Tlmbar_deltaPrime_minus;         // T_l\bar{m}^{(δ',-)}

            std::span<const Complex> rho_at_rmin;                     // rho
            std::span<const Complex> rhoPrimeOverRho_at_rmin;         // rho'/rho
            std::span<const Complex> tau0_at_rmin;                    // tau^\circ
            std::span<const Complex> tau0bar_at_rmin;                 // \bar{tau}^\circ

            std::span<const Real> delta_over_2Sigma_at_rmin;          // Δ/(2Σ)
            std::span<const Real> u_coeff_at_rmin;                    // (r^2+a^2)/Σ
        };

        explicit InitialDataXnn(std::vector<Real> zvals, std::size_t ir_min = 0)
                : zvals_(std::move(zvals)), ir_min_(ir_min)
        {
            if (zvals_.empty()) {
                throw std::runtime_error("InitialDataXnn: zvals cannot be empty.");
            }
        }

        [[nodiscard]] std::size_t size() const noexcept { return zvals_.size(); }
        [[nodiscard]] std::size_t ir_min() const noexcept { return ir_min_; }
        [[nodiscard]] const std::vector<Real>& zvals() const noexcept { return zvals_; }

        [[nodiscard]] Real X0(std::size_t iz, const BoundaryData& bd) const {
            check_sizes_(bd);

            const Complex rho  = bd.rho_at_rmin[iz];
            const Complex rhob = std::conj(rho);

            const Complex rhs =
                    bd.Tln_delta_minus[iz]
                    - bd.delta_over_2Sigma_at_rmin[iz] * bd.Tll_delta_minus[iz]
                    + bd.u_coeff_at_rmin[iz] * bd.du_Tll_deltaPrime_minus[iz]
                    - (rho + rhob) * bd.rhoPrimeOverRho_at_rmin[iz] * bd.Tll_deltaPrime_minus[iz]
                    - (rho  * bd.tedth_Tlm_deltaPrime_minus[iz]
                       + rhob * bd.tedth_Tlmbar_deltaPrime_minus[iz])
                    + Real(2) * rho * rhob *
                      (bd.tau0_at_rmin[iz]    * bd.Tlmbar_deltaPrime_minus[iz]
                       + bd.tau0bar_at_rmin[iz] * bd.Tlm_deltaPrime_minus[iz]);

            const Real denom = (rho + rhob).real();
            if (std::abs(denom) < 1e-14) {
                throw std::runtime_error("InitialDataXnn: denominator (rho+rhob) too small.");
            }

            // x_nn is real; discard tiny imaginary numerical pollution
            return Real(2) * rhs.real() / denom;
        }

        [[nodiscard]] StateVec make_state(std::size_t iz, const BoundaryData& bd) const {
            return StateVec{ X0(iz, bd) };
        }

        [[nodiscard]] std::vector<StateVec> make_all_states(const BoundaryData& bd) const {
            check_sizes_(bd);
            std::vector<StateVec> out;
            out.reserve(zvals_.size());

            for (std::size_t iz = 0; iz < zvals_.size(); ++iz) {
                out.push_back(make_state(iz, bd));
            }
            return out;
        }

        template <typename FieldX>
        void fill_boundary_row(FieldX& Xnn, const BoundaryData& bd) const {
            check_sizes_(bd);
            check_field_shape_(Xnn, "Xnn");

            for (std::size_t iz = 0; iz < zvals_.size(); ++iz) {
                Xnn(ir_min_, iz) = X0(iz, bd);
            }
        }

    private:
        std::vector<Real> zvals_;
        std::size_t ir_min_{0};

        void check_sizes_(const BoundaryData& bd) const {
            const std::size_t Nz = zvals_.size();

            if (bd.Tln_delta_minus.size()               != Nz ||
                bd.Tll_delta_minus.size()               != Nz ||
                bd.du_Tll_deltaPrime_minus.size()       != Nz ||
                bd.Tll_deltaPrime_minus.size()          != Nz ||
                bd.tedth_Tlm_deltaPrime_minus.size()    != Nz ||
                bd.tedth_Tlmbar_deltaPrime_minus.size() != Nz ||
                bd.Tlm_deltaPrime_minus.size()          != Nz ||
                bd.Tlmbar_deltaPrime_minus.size()       != Nz ||
                bd.rho_at_rmin.size()                   != Nz ||
                bd.rhoPrimeOverRho_at_rmin.size()       != Nz ||
                bd.tau0_at_rmin.size()                  != Nz ||
                bd.tau0bar_at_rmin.size()               != Nz ||
                bd.delta_over_2Sigma_at_rmin.size()     != Nz ||
                bd.u_coeff_at_rmin.size()               != Nz)
            {
                throw std::runtime_error(
                        "InitialDataXnn: all boundary arrays must have size zvals.size()."
                );
            }
        }

        template <typename Field>
        void check_field_shape_(const Field& field, const std::string& name) const {
            if (ir_min_ >= field.Nr()) {
                throw std::runtime_error("InitialDataXnn: ir_min out of range for " + name);
            }
            if (field.Nz() != zvals_.size()) {
                throw std::runtime_error("InitialDataXnn: Nz mismatch for " + name);
            }
        }
    };


} // namespace ghz


#endif //GHZ_NUMERIC_INITIALDATA_HPP
