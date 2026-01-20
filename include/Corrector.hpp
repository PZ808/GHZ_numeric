//
// Created by Peter Zimmerman on 31.10.25.
//


#ifndef GHZ_NUMERIC_CORRECTOR_HPP
#define GHZ_NUMERIC_CORRECTOR_HPP

#pragma once

#include "GhzTypes.hpp"
#include "GHPFieldVectorized.hpp"
#include "KerrMetricOutgoing.hpp"
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
     * @class LBuilderBase
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

    class LBuilderBase {
    public:
        // declares a pure virtual function call operator
        virtual StateVec operator()(const StateVec &y, // ref to input state vec
                                    Real r, // radial coord
                                    Real z, // polar coord
                                    const Complex &rho, // kerr rho at this (r,z)
                                    const ghp::HeldBackgroundFieldsVectorized &held, // additional held fields pre-evaluated a z_nodes
                                    size_t iz // index in z direction
        ) const = 0; // all derived classes must implement this

        virtual ~LBuilderBase() = default;
    };


/**
    ** @class LBuilderX_mmbar
     * @brief Implementation of LBuilderBase for the X_mmbar level.
     *
     * @details This class computes the left-hand-side (LHS) operator for the X_mmbar level, \n
     * which means explicitly L[X] in * pd_r^2[X_mmbar] = L[Xmmbar] + T_{ll} transport equation \n
     *  rho^2 pd_r (rhob^2/rho^3 pd_r [rho*X_mmbar/rhob] ) = T_{ll} \n
     *  giving L = (rho+rhob) pd_r + (rho - rhob)^2 \n \n
     * which involves a 4-component real vector layout: [ReX, ReX', ImX, ImX']. \n
     * It defines the behavior of the LHS operator for this specific level. \n
     */
    class LBuilderXmmbar : public LBuilderBase {

    public:
        StateVec operator()(const StateVec &y,
                            Real r, Real z, const Complex &rho,
                            const ghp::HeldBackgroundFieldsVectorized &held,
                            size_t) const override {
            StateVec dy(4, 0.0);
            Complex rhob = std::conj(rho);

            dy[0] = y[1]; // d(ReX)/dr
            dy[2] = y[3]; // d(ImX)/dr

            dy[1] = ((rho + rhob) * y[1] + (rho - rhob) * (rho - rhob) * y[0]).real();
            dy[3] = ((rho + rhob) * y[3] + (rho - rhob) * (rho - rhob) * y[2]).imag(); // should be zero for X_mmbar

            return dy;
        }
    };
    /**
       * @class LBuilderXnm
       * @brief Concrete implementation of LBuilderBase for the X_nm level. \n
       *
       * This class computes the LHS operator for the X_nm level, for the eqn.  \n
       * rho/(2(rho + rhob)) * d/dr ( (rho + rhob)^2 * d/dr [X_nm/(rho*(rho + rhob))] ) = T_lm + N[x_mmbar] \n
       * L[f] = rho dX/dr + rho^2 f - rhob*(rho-rhob)f
       * where N[x_mmbar] is a precomputed __linear__ operator applied to x_mmbar. \n
       */
    class LBuilderXnm : public LBuilderBase {

    public :
        LBuilderXnm(const GHPSpectral& NX_mmbar,
                    const RVector& r_grid,
                    const RVector& z_grid)
                : F_field_(std::move(NX_mmbar)), // forcing term from X_mmbar
                  r_grid_(r_grid), z_grid_(z_grid) {}

        ode::StateVec operator()(const StateVec &y,
                                 Real r, Real /*z*/,
                                 const Complex &rho,
                                 const ghp::HeldBackgroundFieldsVectorized & /*held*/,
                                 size_t iz) const override {
            // Expect y.size() == 4 (ReX, ReX', ImX, ImX')
            ode::StateVec dy(4, 0.0);

            // unpack X_nm into complex quantities
            Complex Xnm(y[0], y[2]);
            Complex Vnm(y[1], y[3]);
            Real z = z_grid_[iz];

            Complex rhob = std::conj(rho);

            Complex LXnm = 2.0_r * (rho * rho * Xnm - rhob * (rho - rhob) * Xnm + rho * Vnm)
                           + 2.0_r * F_field_(findNearestRIndex(r, r_grid_), iz).value(); // \pd_r^2 X_{nm}

            // pack dy
            dy[0] = std::real(Vnm);
            dy[2] = std::imag(Vnm);
            dy[1] = std::real(LXnm);
            dy[3] = std::imag(LXnm);

            return dy;
        }

    private:

        const GHPSpectral F_field_;
        const RVector r_grid_{};
        const RVector z_grid_{};

    };// class LBuilderXnm



    /**
      * @class LBuilderXnn
      * @brief Concrete implementation of LBuilderBase for X_nn  \n
      *
      * This class computes the LHS operator for the X_nn eqn.  \n
      *  1/2 * (rho+rhob)^2 * Thorn[ X_nn/(rho+rhob) ] = Tln + Re U[X_mmbar] + Re V[X_nm] \n
      * L[f] = - rho*rhob/r^2 * ( r^2-a^2*z^2 ) * f -  Sigma/r( Re U[X_mmbar] + Re V[X_nm] )
      * where U[x_mmbar] is a precomputed __linear__ operator applied to x_mmbar \n
      * and V[x_nm] is a precomputed __linear__ operator applied to x_nm. \n
      */
    class LBuilderXnn : public LBuilderBase {
    public:
        LBuilderXnn(const GHPSpectral &UX_mmbar, const GHPSpectral &VX_nm,
                    const RVector &r_grid, const RVector &z_grid)
                : UX_field_(UX_mmbar),
                  VX_field_(VX_nm),
                  r_grid_(r_grid),
                  z_grid_(z_grid) {}

        ode::StateVec operator()(const ode::StateVec &y,
                                 Real r, Real /*z*/,
                                 const Complex &rho,
                                 const ghp::HeldBackgroundFieldsVectorized & /*held*/,
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

    }; // class LBuilderXnn

/**
 * @class ZSliceSolver
 * @brief Integrates the GHZ transport ODE system along the radial grid for a fixed z-slice.
 *
 * This class performs the radial solve: \n
 *
 *      dy/dr = L(y, r, z_fixed) + S(r, z_fixed) \n
 *
 * where:
 *   - L is the left-hand-side operator constructed externally (e.g. from LBuilder_X) \n
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
                     const ghp::HeldBackgroundFieldsVectorized &held,
                     bool use_chebyshev = false)
                : metric_(metric), r_grid_(r_grid), Nz_(Nz), held_(held), use_chebyshev_(use_chebyshev) {
            buildZGrid();
        }

        void setSourceData(const std::vector<std::vector<Real>> &source) { source_data_ = source; }

        void setLBuilder(const LBuilderBase &builder) { L_builder_ptr_ = &builder; }

        StateMat solve()
        {
            StateMat sol_2D(Nz_);

#pragma omp parallel for
            for (size_t iz = 0; iz < Nz_; ++iz) {
                sol_2D[iz] = solveSingleZ(iz);
            }
            return sol_2D;
        }

    private:
        void buildZGrid() {
            if (use_chebyshev_)
                z_grid_ = spectral::SpectralDiffer::build_chebyshev_lobatto_nodes(static_cast<int>(Nz_));
            else
                z_grid_ = spectral::SpectralDiffer::build_legendre_gauss_lobatto_nodes(static_cast<int>(Nz_));
        }

        std::vector<ode::StateVec> solveSingleZ(size_t iz) {
            Real z = z_grid_[iz];

            ode::LHSFunc lhs = [this, iz](const ode::StateVec &y, Real r, Real z_unused) {
                Complex rho = -Real(1) / (r - Complex(0.0, metric_.a() * z_grid_[iz]));
                return (*L_builder_ptr_)(y, r, z_unused, rho, held_, iz);
            };

            ode::SourceFunc src = [this, iz](Real, Real) {
                ode::StateVec S(Nfields_, 0.0);
                for (int i = 0; i < Nfields_; ++i) S[i] = source_data_[iz][i];
                return S;
            };

            ode::TransportODESystem sys{lhs, src};
            return ode::solve_single_z(r_grid_, y0_, sys, z, abs_tol_, rel_tol_, initial_step_);
        }

        const KerrMetricOutgoing &metric_;
        const std::vector<Real> &r_grid_;
        size_t Nz_;
        bool use_chebyshev_= false;
        const ghp::HeldBackgroundFieldsVectorized &held_;
        std::vector<Real> z_grid_;
        const LBuilderBase *L_builder_ptr_ = nullptr;
        std::vector<std::vector<Real>> source_data_;
        int Nfields_ = 4;
        ode::StateVec y0_;
        Real abs_tol_ = 1e-12;
        Real rel_tol_ = 1e-12;
        Real initial_step_ = 1e-4;
    }; // class ZSliceSolver




    class SpectralFieldBuilder {
    public:

        static std::map<std::string, spectral::SpectralGHPVectorized>//std::vector<spectral::SpectralGHPVectorized>
        buildFromSolution(const std::vector<std::vector<ode::StateVec>> &sol_2D,
                          const std::vector<Real> &r_grid,
                          const std::vector<Real> &z_grid,
                          int p, int q, int m, int k, int n, Real a) {
            // declare map of fields to build
            std::map <std::string, spectral::SpectralGHPVectorized> field_map;

            GHPSpectral X(r_grid.size(), z_grid.size(), {m, k, n},
                          GHPScalar<Complex>(teuk::zeroC, p, q));
            GHPSpectral thorn_X(r_grid.size(), z_grid.size(), {m, k, n},
                                      GHPScalar<Complex>(teuk::zeroC, p + 1, q + 1));
            GHPSpectral edth_X(r_grid.size(), z_grid.size(), {m, k, n},
                               GHPScalar<Complex>(teuk::zeroC, p - 1, q + 1));
            GHPSpectral edthBar_X(r_grid.size(), z_grid.size(), {m,k,n},
                                  GHPScalar<Complex>(teuk::zeroC, p + 1, q - 1));

#pragma omp parallel for collapse(2) default(none) shared(X, sol_2D, z_grid, r_grid, p, q)
            for (size_t iz = 0; iz < z_grid.size(); ++iz) {
                for (size_t ir = 0; ir < r_grid.size(); ++ir) {
                    Complex val(sol_2D[iz][ir][0], sol_2D[iz][ir][2]);
                    X.set_index(ir, iz, GHPScalar<Complex>(val, p, q));
                }
            }

            // edth / edthBar computation
#pragma omp parallel for default(none) shared(X, edth_X, edthBar_X, r_grid, z_grid, p, q, a)
            for (size_t ir = 0; ir < r_grid.size(); ++ir) {
                spectral::SpectralDiffer differ(z_grid.size(), r_grid.size());
                auto slice = X.slice_R(ir);
                auto slice_edthH = edth_X.slice_R(ir);
                auto slice_edthBarH = edthBar_X.slice_R(ir);

                differ.edthH_inplace_RSliceV(slice, slice_edthH, a);
                differ.edthBarH_inplace_RSliceV(slice, slice_edthBarH, a);
                for (size_t iz = 0; iz < z_grid.size(); ++iz) {
                    edth_X.set_index(ir, iz, GHPScalar<Complex>(slice_edthH[iz].value(), p - 1, q + 1));
                    edthBar_X.set_index(ir, iz, GHPScalar<Complex>(slice_edthBarH[iz].value(), p + 1, q - 1));
                }
            }

            field_map.emplace("X", X);
            field_map.emplace("thX", thorn_X);
            field_map.emplace("edthX", edth_X);
            field_map.emplace("edthBarX", edthBar_X);

            return field_map;
        }
    };

/**
 * @class KerrODECorrector
 * @brief Solves hierarchical GHZ transport ODEs along r for each fixed-z slice.
 *
 * This class provides a reusable driver for solving radial ODE systems of the form
 *
 *      (D + N) X(r, z) = S(r, z)
 *
 * for each Legendre/Chebyshev z-node. It is designed for the GHZ transport hierarchy
 * where one solves:
 *
 *   1. X_mmbar   (lowest level)
 *   2. X_nm      (uses derivatives of X_mmbar as source)
 *   3. X_nn      (uses derivatives of X_nm as source)
 *
 * At each z index (iz), the class:
 *   - evaluates background fields (rho, τ°, etc.) using KerrMetricOutgoing + Held fields
 *   - constructs the left-hand-side operator using the user-supplied L_builder functor
 *   - injects the appropriate source vector from setSourceData()
 *   - calls ode::solve_single_z to integrate across the radial grid r_grid_
 *
 * The class works with Chebyshev-Lobatto or Legendre-Gauss-Lobatto z-grids.
 *
 * Typical usage:
 *  - Construct KerrODECorrector(metric, r_grid, Nz, held)
 *  - Provide the radial source values for each z using setSourceData()
 *  - Provide the operator callback using setLBuilder()
 *  - Call solve() to obtain a 2D array sol[z][ir] of ode::StateVec
 *
 * The solver does NOT own the metric or Held field objects; they must outlive it.
 *
 * ### Constructor Arguments
 * @param metric    Reference to KerrMetricOutgoing (defines Δ, ρ, spin a, etc.).
 * @param r_grid    Radial nodes (Chebyshev or user-supplied).
 * @param Nz        Number of z collocation points.
 * @param held      HeldBackgroundFieldsVectorized providing τ°, ρ°, and other NP fields.
 * @param use_chebyshev If true, z-grid uses Chebyshev-Lobatto; otherwise LGL.
 *
 * ### Required User Callbacks
 * setLBuilder(): must supply a function constructing LHS operator
 *      ode::StateVec L(y, r, z, rho, held_fields, iz)
 *
 * setSourceData(): supply S(z, i) for each field component.
 *
 * ### Example
 *
 * @code
 * // 1. Construct solver
 * KerrODECorrector corr(metric, r_grid, Nz, held_background);
 *
 * // 2. Fill source data for all z slices (3 fields here)
 * std::vector<std::vector<Real>> src(Nz, std::vector<Real>(3));
 * for (size_t iz = 0; iz < Nz; ++iz) {
 *     src[iz][0] = S_mmbar(iz);   // example
 *     src[iz][1] = S_nm(iz);
 *     src[iz][2] = S_nn(iz);
 * }
 * corr.setSourceData(src);
 *
 * // 3. Provide LHS operator (GHZ N-operator for each sector)
 * corr.setLBuilder(
 *     [&](const ode::StateVec& X, Real r, Real z, const Complex& rho,
 *         const ghp::HeldBackgroundFieldsVectorized& held, const size_t& iz)
 *     {
 *         ode::StateVec out(3);
 *         // Example dummy structure:
 *         out[0] = compute_L_mmbar(X, r, z, rho, held, iz);
 *         out[1] = compute_L_nm   (X, r, z, rho, held, iz);
 *         out[2] = compute_L_nn   (X, r, z, rho, held, iz);
 *         return out;
 *     }
 * );
 *
 * // 4. Solve full 2D (z × r) problem
 * auto solution = corr.solve();
 *
 * // solution[iz][ir][field_index] contains result for each component.
 * @endcode
 */
    class KerrODECorrector {
    public:
        KerrODECorrector(const KerrMetricOutgoing &metric,
                         const std::vector<Real> &r_grid,
                         size_t Nz,
                         const ghp::HeldBackgroundFieldsVectorized &held,
                         bool use_chebyshev = true)
                : metric_(metric), r_grid_(r_grid), Nz_(Nz), held_(held), use_chebyshev_(use_chebyshev) {
            buildZGrid();
            Nfields_ = 3; // or adapt based on source_data
            y0_ = ode::StateVec(Nfields_, 0.0);
        }

        void setSourceData(const std::vector<std::vector<Real>> &source) { source_data_ = source; }

        void setLBuilder(const std::function<ode::StateVec(const ode::StateVec &, Real, Real, const Complex &,
                                                           const ghp::HeldBackgroundFieldsVectorized &,
                                                           const size_t &)> &L) {
            L_builder_ = L;
        }

        // Solve along all z-slices
        std::vector<std::vector<ode::StateVec>> solve() {
            std::vector<std::vector<ode::StateVec>> sol_2D;
            sol_2D.reserve(Nz_);

            for (size_t iz = 0; iz < Nz_; ++iz) {
                sol_2D.push_back(solveSingleZ(iz));
            }

            return sol_2D;
        }

    private:
        std::vector<ode::StateVec> solveSingleZ(size_t iz) {
            Real z = z_grid_[iz];

            ode::LHSFunc lhs = [this, iz](const ode::StateVec &y, Real r, Real) {
                Complex rho = -Real(1) / (r - Complex(0.0, metric_.a() * z_grid_[iz]));
                return L_builder_(y, r, z_grid_[iz], rho, held_, iz);
            };

            ode::SourceFunc src = [this, iz](Real, Real) {
                ode::StateVec S(Nfields_);
                for (int i = 0; i < Nfields_; ++i) {
                    S[i] = source_data_[iz][i];
                }
                return S;
            };

            ode::TransportODESystem sys{lhs, src};
            return ode::solve_single_z(r_grid_, y0_, sys, z, abs_tol_, rel_tol_, initial_step_);
        }

        void buildZGrid() {
            if (use_chebyshev_) {
                z_grid_ = spectral::SpectralDiffer::build_chebyshev_lobatto_nodes(static_cast<int>(Nz_));
            } else {
                z_grid_ = spectral::SpectralDiffer::build_legendre_gauss_lobatto_nodes(static_cast<int>(Nz_));
            }
        }

    private:
        const KerrMetricOutgoing &metric_;
        const std::vector<Real> &r_grid_;
        size_t Nz_;
        int Nfields_;
        bool use_chebyshev_;
        ode::StateVec y0_;
        std::vector<Real> z_grid_;
        std::vector<std::vector<Real>> source_data_;
        const ghp::HeldBackgroundFieldsVectorized &held_;
        std::function<ode::StateVec(const ode::StateVec &, Real, Real, const Complex &,
                                    const ghp::HeldBackgroundFieldsVectorized &, const size_t &)> L_builder_;
        Real abs_tol_ = 1e-12;
        Real rel_tol_ = 1e-12;
        Real initial_step_ = 1e-4;
    };


/**
 * @class GHZTransportSolver
 * @brief High-level driver for solving complex GHZ transport ODEs for each z-slice.
 *
 * This class replaces the old procedural routine `corrector_solve(...)` with an
 * object-oriented interface. It solves hierarchical GHZ transport equations along
 * the radial grid r for every fixed z collocation point.
 *
 * The solver performs, for each z-slice:
 *   - evaluates the z-grid (Chebyshev or LGL)
 *   - computes ρ(r, z) = -1 / (r - i a z)
 *   - constructs the left-hand-side operator via the user-supplied L_builder
 *   - injects the source vector source_data[iz][field]
 *   - calls ode::solve_single_z(...) to integrate the radial ODE system
 *   - stores a 2D solution array sol[z][r]
 *
 * This is used for all three GHZ transport levels:
 *     X_mmbar  →  X_nm  →  X_nn
 *
 * ### Constructor arguments
 * @param kerr        Reference to KerrMetricOutgoing (defines geometry, spin a).
 * @param r_grid      Radial collocation grid.
 * @param Nz          Number of z collocation points.
 * @param held        Vectorized Held background fields (ρ°, τ°, etc.).
 * @param abs_tol     Absolute tolerance for ODE integrator.
 * @param rel_tol     Relative tolerance for ODE integrator.
 * @param initial_step Initial step size for the integrator.
 * @param use_cheb_z  If true: Chebyshev-Lobatto z-grid; otherwise LGL.
 *
 * ### Usage example
 * @code
 * GHZTransportSolver solver(
 *      metric, r_grid, Nz, held_fields,
 *      1e-12, 1e-12, 1e-4,
 *      true   // use Chebyshev in z
 * );
 *
 * // LHS operator builder
 * GHZTransportSolver::LBuilderFunc Lb =
 *     [&](const ode::StateVec& y, Real r, Real z, const Complex rho,
 *         const ghp::HeldBackgroundFieldsVectorized& held, size_t iz)
 *     {
 *         ode::StateVec out(y.size());
 *         out[0] = L_mmbar(y, r, z, rho, held, iz);
 *         out[1] = L_nm(y, r, z, rho, held, iz);
 *         out[2] = L_nn(y, r, z, rho, held, iz);
 *         return out;
 *     };
 *
 * // Source for each z and each GHZ component
 * std::vector<std::vector<Real>> source(Nz, std::vector<Real>(3));
 * fill_sources_for_level(source);
 *
 * auto sol2D = solver.solve_single_level(Lb, source);
 * // sol2D[iz][ir][component] gives result X_mmbar, X_nm or X_nn
 * @endcode
 *
 */
    class GHZTransportSolver {

    public:
        using LBuilderFunc =
                std::function<ode::StateVec(const ode::StateVec &, Real, Real, const Complex,
                                            const ghp::HeldBackgroundFieldsVectorized &, size_t)>;

        GHZTransportSolver(const KerrMetricOutgoing &kerr,
                           const std::vector<Real> &r_grid,
                           size_t Nz,
                           const ghp::HeldBackgroundFieldsVectorized &held,
                           Real abs_tol,
                           Real rel_tol,
                           Real initial_step,
                           bool use_cheb_z)
                : kerr_(kerr),
                  r_grid_(r_grid),
                  Nz_(Nz),
                  held_(held),
                  abs_tol_(abs_tol),
                  rel_tol_(rel_tol),
                  initial_step_(initial_step),
                  use_cheb_z_(use_cheb_z) {
            // Build z-grid (LGL or Chebyshev)
            if (use_cheb_z_) {
                z_grid_ = spectral::SpectralDiffer::build_chebyshev_lobatto_nodes(Nz_);
            } else {
                z_grid_ = spectral::SpectralDiffer::build_legendre_gauss_lobatto_nodes(Nz_);
            }
        } // constructor

        /// ------------------------------------------------------------------------
        /// Solve a single GHZ hierarchy level.
        ///
        /// Equivalent to your procedural function `corrector_solve(...)`.
        ///
        /// Args:
        ///     L_builder:   functor that builds LHS operator L(y,r,z,rho,held,iz)
        ///     source_data: source_data[iz][field]
        ///
        /// Returns:
        ///     sol_2D[z][r] : the ODE solution at each r for each z.
        /// ------------------------------------------------------------------------
        StateMat solve_single_level(const LBuilderFunc &L_builder,
                           const std::vector<std::vector<Real>> &source_data) {
            const Real a = kerr_.a();

            // Number of fields in system
            const int Nfields = (!r_grid_.empty() ? source_data[0].size() : 3);

            // Zero initial condition
            ode::StateVec y0(Nfields, 0.0);

            // Storage for [z][r]
            std::vector<std::vector<ode::StateVec>> sol_2D;
            sol_2D.reserve(Nz_);

            // Loop over z-slices
            for (size_t iz = 0; iz < Nz_; ++iz) {
                Real z = z_grid_[iz];

                // Build LHS function for this z
                ode::LHSFunc lhs = [&](const ode::StateVec &y, Real r, Real z_unused) {
                    Complex rho = -Real(1) / (r - Complex(0.0, a * z));
                    return L_builder(y, r, z_unused, rho, held_, iz);
                };

                // Build source function for this z
                ode::SourceFunc src = [&](Real, Real) {
                    ode::StateVec S(Nfields);
                    for (int i = 0; i < Nfields; ++i)
                        S[i] = source_data[iz][i];
                    return S;
                };

                ode::TransportODESystem sys{lhs, src};

                // Solve along r-grid
                sol_2D.push_back(
                        ode::solve_single_z(r_grid_, y0, sys, z,
                                            abs_tol_, rel_tol_, initial_step_)
                );
            }

            return sol_2D;
        } // solve_single_level

        /// Direct access to the z-grid
        const std::vector<Real> &z_grid() const { return z_grid_; }

    private:
        const KerrMetricOutgoing &kerr_;
        const std::vector<Real> &r_grid_;
        size_t Nz_;
        const ghp::HeldBackgroundFieldsVectorized &held_;

        Real abs_tol_;
        Real rel_tol_;
        Real initial_step_;
        bool use_cheb_z_;

        std::vector<Real> z_grid_;
    }; // class GHZTransportSolver


    // Function to read source data from a file
    // The file is assumed to contain Nfields columns per line for each r-grid point
    std::vector<std::vector<Real>> read_source_file(const std::string &filename, size_t Nr, size_t Nfields);

    // Corrector solver for coupled ODEs
    std::vector<std::vector<ode::StateVec>> corrector_solve(
            const KerrMetricOutgoing &kerr_out,
            const std::vector<Real> &r_grid,
            int Nz,
            const ghp::HeldBackgroundFieldsVectorized &held,
            const std::function<ode::StateVec(const ode::StateVec &, Real, Real)> &LHS_builder,
            const std::vector<std::vector<Real>> &source_data,
            Real abs_tol = 1e-9,
            Real rel_tol = 1e-9,
            Real initial_step = 1e-3,
            bool use_chebyshev = true
    );

} // namespace ghz


#endif //GHZ_NUMERIC_CORRECTOR_HPP