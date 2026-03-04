//
// Created by Peter Zimmerman on 31.10.25.
//


#ifndef GHZ_NUMERIC_CORRECTOR_HPP
#define GHZ_NUMERIC_CORRECTOR_HPP

#pragma once

#include "GhzTypes.hpp"
#include "GHPFieldVectorized.hpp"
#include "KerrMetricOutgoing.hpp"
#include "KinnersleyTetrad.hpp"
#include "KinnersleySpectralHeldOperators.hpp"
#include "ODE.hpp"
#include "HeldScalars.hpp"
#include "SpectralDiffer.hpp"  // your node builders


#include <utility>
#include <vector>
#include <string>

namespace ghz {
    using namespace teuk::literals;
    using StateVec =  ode::StateVec;
    using StateMat =  std::vector<std::vector<StateVec>>;
    using RVector = std::vector<Real>;
    using CVector = std::vector<Complex>;
    using GHPSpectral = spectral::SpectralGHPVectorized;
    /**
     * @class BuilderBase
     * @brief Abstract base class for building LHS operators for different components of the corrector.
     */
    /**
      * @brief Finds the nearest radial grid index for a given radial coordinate.
      * @param r The radial coordinate.
      * @return The index of the nearest radial grid point.
      */
    size_t findNearestRIndex(Real r, const RVector& r_grid) {
        size_t ir = 0;
        Real min_diff = std::abs(r - r_grid[0]);
        for (size_t i = 1; i < r_grid.size(); ++i) {
            Real d = std::abs(r - r_grid[i]);
            if (d < min_diff) {
                min_diff = d;
                ir = i;
            }
        }
        return ir;
    }

    class BuilderBase {
    public:
        // declares a pure virtual function call operator
        virtual StateVec operator()(const StateVec &y, // ref to input state vec
                                    Real r, // radial coord
                                    Real z, // polar coord
                                    const Complex &rho, // kerr rho at this (r,z)
                                    const ghp::HeldBackgroundFieldsVectorized<CoordType> &held, // additional held fields pre-evaluated a z_nodes
                                    size_t iz // index in z direction
        ) const = 0; // all derived classes must implement this

        virtual ~BuilderBase() = default;
    };


/**
    ** @class BuilderX_mmbar
     * @brief Implementation of BuilderBase for the X_mmbar level.
     *
     * @details This class computes the lower-order "L-operator" for the X_mmbar level, \n
     * which means operationally \n
     * keep the highest r-derivative on the left and move everything else to the RHS  \n
     * and divide so that the coefficient of the highest derivative term is one.. \n
     * L[X] in * pd_r^2[X_mmbar] = L[Xmmbar] + T_{ll} transport equation \n
     *  rho^2 pd_r (rhob^2/rho^3 pd_r [rho*X_mmbar/rhob] ) = T_{ll} \n
     *  giving L = (rho+rhob) pd_r + (rho - rhob)^2 \n
     * which involves a 4-component real vector layout: [ReX, ReX', ImX, ImX']. \n
     * It defines the behavior of the LHS operator for this specific level. \n
     */
    class BuilderXmmbar : public BuilderBase {

    public:
        StateVec operator()(const StateVec &y,
                            Real r, Real z, const Complex &rho,
                            const ghp::HeldBackgroundFieldsVectorized<CoordType> &held,
                            size_t) const override {

            StateVec dy(4, 0.0); // initially zero. fill with LHS deriv operator applied to y

            Complex rhob = std::conj(rho);

            dy[0] = y[1]; // d(ReX)/dr
            dy[2] = y[3]; // d(ImX)/dr

            dy[1] = (rho + rhob).real() * y[1] + ((rho - rhob) * (rho - rhob)).real() * y[0]; // d^2 ReX/dr^2 = L[X_mmbar] applied to ReX
            dy[3] = (rho + rhob).real() * y[3] + ((rho - rhob) * (rho - rhob)).real() * y[2]; // d^2 ImX/dr^2 = L[X_mmbar] applied to ImX

            return dy;
        }
    };

    /**
       * @class BuilderXnm
       * @brief Concrete implementation of BuilderBase for the X_nm level. \n
       *
       * This class computes the operator for the X_nm level, for the eqn.  \n
       * rho/(2(rho + rhob)) * d/dr ( (rho + rhob)^2 * d/dr [X_nm/(rho*(rho + rhob))] ) = T_lm + N[x_mmbar] \n
       * L[f] = rho dX/dr + rho^2 f - rhob*(rho-rhob)f
       * where N[x_mmbar] is a precomputed __linear__ operator applied to x_mmbar. \n
       */
    class BuilderXnm : public BuilderBase {

    public :
        BuilderXnm(const GHPSpectral& NX_mmbar,
                    const RVector& r_grid,
                    const RVector& z_grid)
                : NX_source_(std::move(NX_mmbar)), // forcing term from X_mmbar
                  r_grid_(r_grid), z_grid_(z_grid) {}

        ode::StateVec operator()(const StateVec &y,
                                 Real r, Real /*z*/,
                                 const Complex &rho,
                                 const ghp::HeldBackgroundFieldsVectorized<CoordType> & /*held*/,
                                 size_t iz) const override {
            // Expect y.size() == 4 (ReX, ReX', ImX, ImX')
            ode::StateVec dy(4, 0.0);

            // unpack X_nm into complex quantities
            Complex Xnm(y[0], y[2]); // ReXnm and ImXnm state vector
            Complex Vnm(y[1], y[3]); // dReX_nm/dr and dImXnm/dr
            Real z = z_grid_[iz];

            Complex rhob = std::conj(rho);

            Complex LXnm = 2.0_r * (rho * rho * Xnm - rhob * (rho - rhob) * Xnm + rho * Vnm)
                           + 2.0_r * NX_source_(findNearestRIndex(r, r_grid_), iz).value(); // \pd_r^2 X_{nm}

            // pack dy
            dy[0] = std::real(Vnm);
            dy[2] = std::imag(Vnm);
            dy[1] = std::real(LXnm);
            dy[3] = std::imag(LXnm);

            return dy;
        }
    private:

        const GHPSpectral NX_source_;
        const RVector r_grid_{};
        const RVector z_grid_{};

    };// class BuilderXnm



    /**
      * @class BuilderXnn
      * @brief Concrete implementation of BuilderBase for X_nn  \n
      *
      * This class computes the LHS operator for the X_nn eqn.  \n
      *  1/2 * (rho+rhob)^2 * Thorn[ X_nn/(rho+rhob) ] = Tln + Re U[X_mmbar] + Re V[X_nm] \n
      * L[f] = - rho*rhob/r^2 * ( r^2-a^2*z^2 ) * f -  Sigma/r( Re U[X_mmbar] + Re V[X_nm] )
      * where U[x_mmbar] is a precomputed __linear__ operator applied to x_mmbar \n
      * and V[x_nm] is a precomputed __linear__ operator applied to x_nm. \n
      */
    class BuilderXnn : public BuilderBase {
    public:
        BuilderXnn(const GHPSpectral &UX_mmbar, const GHPSpectral &VX_nm,
                    const RVector &r_grid, const RVector &z_grid)
                : UX_field_(UX_mmbar),
                  VX_field_(VX_nm),
                  r_grid_(r_grid),
                  z_grid_(z_grid) {}

        ode::StateVec operator()(const ode::StateVec &y,
                                 Real r, Real /*z*/,
                                 const Complex &rho,
                                 const ghp::HeldBackgroundFieldsVectorized<CoordType> & /*held*/,
                                 size_t iz) const override {

            // State: y = {Re Xnn, Im Xnn}
            Complex Xnn(y[0], y[1]);

            Complex rhob = std::conj(rho);
            Real z = z_grid_[iz];

            // source term S = T_ln + Re U + Re V
            Complex U = UX_field_(findNearestRIndex(r, r_grid_), iz).value();
            Complex V = VX_field_(findNearestRIndex(r, r_grid_), iz).value();
            Real S = std::real(U + V);   // assuming T_ln already included into U or V

            // compute dXnn/dr from first-order formula
            Complex dX =
                    Complex(2.0 * S, 0.0) / (rho + rhob)
                    + ((rho * rho + rhob * rhob) / (rho + rhob)) * Xnn;

            ode::StateVec dy(2);
            dy[0] = std::real(dX);
            dy[1] = std::imag(dX);
            return dy;
        }

    private:
        const spectral::SpectralGHPVectorized UX_field_, VX_field_;
        const std::vector<Real> r_grid_, z_grid_;

    }; // class BuilderXnn

/**
 * @class ZSliceSolver
 * @brief Integrates the GHZ transport ODE system along the radial grid for a fixed z-slice.
 *
 * This class performs the radial solve: \n
 *
 *      dy/dr = L(y, r, z_fixed) + S(r, z_fixed) \n
 *
 * where:
 *   - L is the left-hand-side operator constructed externally (e.g. from Builder_X) \n
 *   - S is the vector-valued source term evaluated at the current radial point \n
 *   - y0 is the initial condition at r_min \n
 *
 * The solver operates only along r (1D ODE), and is used by higher-level classes \n
 * such as GHZTransportSolver to build the full 2D (z × r) solution. \n
 *
 * Recommended solver tolerances for GHZ transport:
 *   abs_tol ≈ 1e-12
 *   rel_tol ≈ 1e-12
 */

    class ZSliceSolver {
    public:
        ZSliceSolver(const KerrMetricOutgoing &metric,
                     const RVector &r_grid,
                     size_t Nz,
                     const ghp::HeldBackgroundFieldsVectorized<CoordType> &held,
                     bool use_chebyshev = false)
                : metric_(metric), r_grid_(r_grid), Nz_(Nz), held_(held), use_chebyshev_(use_chebyshev) {
            buildZGrid();
        }

        StateMat solve()
        {
            StateMat sol_2D(Nz_);

#pragma omp parallel for
            for (size_t iz = 0; iz < Nz_; ++iz) {
                sol_2D[iz] = solveSingleZ(iz);
            }
            return sol_2D;
        }

        void setInitialCondition(const ode::StateVec& y0) {
            y0_ = y0;
            Nfields_ = static_cast<int>(y0_.size());
        }

        void setSourceData(const std::vector<std::vector<Real>> &source) { source_data_ = source; }
        void setSourceFunc(ode::SourceFunc f) {
            source_func_ = std::move(f);
        }

        void setBuilder(const BuilderBase &builder) { builder_ptr_ = &builder; }

    private:
        void buildZGrid() {
            if (use_chebyshev_)
                z_grid_ = spectral::SpectralDiffer::build_chebyshev_lobatto_nodes(static_cast<int>(Nz_));
            else
                z_grid_ = spectral::SpectralDiffer::build_legendre_gauss_lobatto_nodes(static_cast<int>(Nz_));
        }
        std::vector<ode::StateVec> solveSingleZ(size_t iz) {
            ode::LHSFunc lhs = [this, iz](const ode::StateVec &y, Real r, Real z_unused) {
                Complex rho = -Real(1) / (r - Complex(0.0, metric_.a() * z_grid_[iz]));
                return (*builder_ptr_)(y, r, z_unused, rho, held_, iz);
            };

            ode::SourceFunc src = [this, iz](Real r, Real /*z_unused*/) {
                if (!source_func_) return ode::StateVec(Nfields_, 0.0);
                return source_func_(r, z_grid_[iz]); // <-- depends on r and z
            };

            ode::TransportODESystem sys{lhs, src};
            return ode::solve_single_z(r_grid_, y0_, sys, z_grid_[iz], abs_tol_, rel_tol_, initial_step_);
        }


        const KerrMetricOutgoing &metric_;
        const std::vector<Real> &r_grid_;
        size_t Nz_;
        bool use_chebyshev_= false;
        const ghp::HeldBackgroundFieldsVectorized<CoordType> &held_;
        std::vector<Real> z_grid_;
        const BuilderBase *builder_ptr_ = nullptr;
        ode::SourceFunc source_func_ = nullptr;
        std::vector<std::vector<Real>> source_data_;
        int Nfields_ = 4;
        ode::StateVec y0_;
        Real abs_tol_ = 1e-12;
        Real rel_tol_ = 1e-12;
        Real initial_step_ = 1e-4;
    }; // class ZSliceSolver

    class GHPDifferX {
    public:
        // Layout for the 4-real state [ReX, ReX', ImX, ImX']
        struct Layout4 {
            int re = 0;
            int reR = 1;
            int im  = 2;
            int imR = 3;
        };

        static std::map<std::string, spectral::SpectralGHPVectorized>
        diff_solution(const std::vector<std::vector<ode::StateVec>> &state_2D,
                          const std::vector<Real> &r_grid,
                          const std::vector<Real> &z_grid,
                          const ghp::HeldBackgroundFieldsVectorized<OutgoingCoords> &held,
                          const KerrMetricOutgoing &metric,
                          const KinnersleyTetrad<OutgoingCoords> &tetrad,
                          int p, int q,
                          int m, int kr, int kz,
                          Layout4 L) {


            std::map<std::string, GHPSpectral> field_map;
            Real a = metric.a();

            const size_t Nr = r_grid.size();
            const size_t Nz = z_grid.size();

            // Base field X has type (p,q)
            GHPSpectral X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p, q));

            // thorn X (if thorn = ∂r) has type (p+1,q+1)
            GHPSpectral thorn_X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p+1, q+1));

            // thorn' X (if thorn' = ∂t) has type (p-1,q-1)
            GHPSpectral thornPH_X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-1, q-1));

            // Held-edth outputs: make sure these output weights match what your differ actually implements
            // (below assumes edthH : (p,q) -> (p,q-2) and edthBarH : (p,q) -> (p-2,q)
            GHPSpectral EdthH_X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p, q-2));
            GHPSpectral EdthBarH_X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-2, q));

            GHPSpectral EdthH_Thorn_X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p+1, q-1));
            GHPSpectral EdthBarH_Thorn_X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, p-1, q+1));

            GHPSpectral EdthH_EdthBarH_X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, -2, -2));
            GHPSpectral EdthBarH_EdthH_X(Nr, Nz, {m, kr, kz}, GHPScalar<Complex>(teuk::zeroC, -2, -2));


            // Fill X and thorn_X directly from the stored ODE state
#pragma omp parallel for collapse(2) default(none) shared(X, thorn_X, state_2D, Nr, Nz, p, q, L)
            for (size_t iz = 0; iz < Nz; ++iz) {
                for (size_t ir = 0; ir < Nr; ++ir) {
                    const auto &y = state_2D[iz][ir];

                    Complex Xval(  y[L.re],  y[L.im]  );
                    Complex dXdr(  y[L.reR], y[L.imR] );   // ← use previous level’s stored radial derivative

                    X.set_index(ir, iz, GHPScalar<Complex>(Xval, p, q));
                    thorn_X.set_index(ir, iz, GHPScalar<Complex>(dXdr, p+1, q+1));
                }
            } // populate X and thorn_X from ODE state

            // Compute Held derivatives on each r-slice
#pragma omp parallel default(none) shared(Nr, Nz, X, thorn_X, held, tetrad,\
EdthH_X, EdthBarH_X, EdthH_Thorn_X,EdthBarH_Thorn_X,EdthBarH_EdthH_X, \
EdthH_EdthBarH_X, r_grid, z_grid, a, p, q)
            {

                spectral::SpectralDiffer differ(static_cast<int>(Nz), static_cast<int>(Nr));
                KinnersleyHeldOperators<OutgoingCoords> held_ops(differ, held, tetrad);

#pragma omp for
                for (size_t ir = 0; ir < Nr; ++ir) {

                    auto slice_X = X.slice_R(ir);
                    auto slice_thornPH_X = thorn_X.slice_R(ir);
                    auto slice_thorn_X = thorn_X.slice_R(ir);
                    auto slice_EdthH_X = EdthH_X.slice_R(ir);
                    auto slice_EdthHBar_X = EdthBarH_X.slice_R(ir);
                    auto slice_EdthH_Thorn_X = EdthH_Thorn_X.slice_R(ir);
                    auto slice_EdthBarH_Thorn_X = EdthBarH_Thorn_X.slice_R(ir);
                    auto slice_EdthH_EdthBarH_X = EdthH_EdthBarH_X.slice_R(ir);
                    auto slice_EdthBarH_EdthH_X = EdthBarH_EdthH_X.slice_R(ir);

                    held_ops.thornPH_inplace_RSliceV(slice_X, slice_thornPH_X);
                    held_ops.edthH_inplace_RSliceV(slice_X, slice_EdthH_X);
                    held_ops.edthBarH_inplace_RSliceV(slice_X, slice_EdthHBar_X);
                    held_ops.edthH_inplace_RSliceV(slice_EdthHBar_X, slice_EdthH_EdthBarH_X);
                    held_ops.edthBarH_inplace_RSliceV(slice_EdthH_X, slice_EdthBarH_EdthH_X);

                    held_ops.edthH_inplace_RSliceV(slice_thorn_X, slice_EdthH_Thorn_X);
                    held_ops.edthBarH_inplace_RSliceV(slice_thorn_X, slice_EdthBarH_Thorn_X);

                    assert(slice_X[0].p() == p && slice_X[0].q() == q);
                    assert(slice_EdthHBar_X[0].p() == p-2 && slice_EdthHBar_X[0].q() == q);
                    assert(slice_EdthH_X[0].p() == p && slice_EdthH_X[0].q() == q-2);
                    }
                } // loop over r-slices to compute Held derivatives

            field_map.emplace("X", X);
            field_map.emplace("thorn_X", thorn_X);
            field_map.emplace("edth_X", EdthH_X);
            field_map.emplace("edthBar_X", EdthBarH_X);
            field_map.emplace("edth_Thorn_X", EdthH_Thorn_X);
            field_map.emplace("edthBar_Thorn_X", EdthBarH_Thorn_X);
            field_map.emplace("edthH_edthBarH_X", EdthH_EdthBarH_X);
            field_map.emplace("edthBarH_edthH_X", EdthBarH_EdthH_X);

            return field_map;

        } // diff_solution

    }; // class GHPDifferX

} // namespace ghz


#endif //GHZ_NUMERIC_CORRECTOR_HPP//

