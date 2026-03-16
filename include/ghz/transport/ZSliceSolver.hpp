//
// Created by Peter Zimmerman on 05.03.26.
//

#ifndef HS1DATA_DATA_ZSLICESOLVER_HPP
#define HS1DATA_DATA_ZSLICESOLVER_HPP

#pragma once

#include "ghz/transport/Corrector.hpp"
#include "ghz/transport/Builders.hpp"

#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ghz {

    using StateVec = ode::StateVec;
    using StateMat = std::vector<std::vector<StateVec>>;
    using RVector  = std::vector<Real>;
/**
 * @class ZSliceSolver
 * @brief Solves a transport ODE system independently on each fixed-z slice of the 2D spectral grid.
 *
 * @details
 * This class performs the radial integration step of the GHZ transport hierarchy. \n
 * For each angular (collocation point) z = z_iz, it solves an ODE system of the form\n
 *\n
 *   dy/dr = L(y; r, z_iz) + S(r, z_iz), \n
 *\n
 * where:\n
 *   - L is the homogeneous left-hand-side operator supplied by a concrete BuilderBase
 *     implementation,\n
 *   - S is the source term supplied through a SourceFunc,\n
 *   - y(r_min, z_iz) = y0(z_iz) is the boundary data at the inner radial boundary.\n
 *\n
 * The key point is that the initial data generally depends on the z-slice, so the solver\n
 * accepts a function initial_data_func(iz) which provides the correct state vector for \n
 * each z index separately.\n
 *\n
 * Operationally, the class:\n
 *   1. loops over all z slices,\n
 *   2. builds the slice-local LHS and source functions,\n
 *   3. integrates along the full radial grid r_grid_,\n
 *   4. returns the full collection of radial solutions as a z-indexed matrix of states.\n
 *
 * This class is agnostic about which corrector level is being solved. It can be used for:\n
 *   - X_mmbar  (real, second-order system),\n
 *
 *   - X_nm     (complex, second-order system written in real form),\n
 *   - X_nn     (real, first-order system),\n
 *   \n
 * provided that the supplied BuilderBase, source function, and per-slice initial data\n
 * all use consistent state-vector conventions.\n
 *
 * The returned object has layout\n
 *\n
 *
 *   sol_2D[iz][ir] = state vector at (r_ir, z_iz),\n
 *
 * so each z-slice stores the full radial evolution of its transport state. \n
 */
    class ZSliceSolver {
    public:
        using InitialDataFunc = std::function<StateVec(size_t)>;

        ZSliceSolver(const KerrMetricOutgoing& metric,
                     const RVector& r_grid,
                     size_t Nz,
                     const ghp::HeldBackgroundFieldsVectorized<CoordType>& held,
                     bool use_chebyshev = false)
                : metric_(metric),
                  r_grid_(r_grid),
                  Nz_(Nz),
                  held_(held),
                  use_chebyshev_(use_chebyshev)
        {
            if (r_grid_.empty()) {
                throw std::invalid_argument("ZSliceSolver: r_grid is empty.");
            }
            if (Nz_ == 0) {
                throw std::invalid_argument("ZSliceSolver: Nz must be > 0.");
            }
            buildZGrid();
        }

        void setBuilder(const BuilderBase& builder) {
            builder_ptr_ = &builder;
        }

        void setSourceFunc(ode::SourceFunc f) {
            source_func_ = std::move(f);
        }

        [[nodiscard]] const RVector& z_grid() const noexcept {
            return z_grid_;
        }

        [[nodiscard]] const RVector& r_grid() const noexcept {
            return r_grid_;
        }

        StateMat solve(const InitialDataFunc& initial_data_func) const {
            if (!builder_ptr_) {
                throw std::runtime_error("ZSliceSolver::solve: builder not set.");
            }
            if (!initial_data_func) {
                throw std::runtime_error("ZSliceSolver::solve: initial_data_func not set.");
            }

            StateMat sol_2D(Nz_);

#pragma omp parallel for
            for (size_t iz = 0; iz < Nz_; ++iz) {
                const StateVec y0 = initial_data_func(iz);
                sol_2D[iz] = solveSingleZ(iz, y0);
            }

            return sol_2D;
        }

    private:
        void buildZGrid() {
            if (use_chebyshev_) {
                z_grid_ = spectral::SpectralDiffer::build_chebyshev_lobatto_nodes(
                        static_cast<int>(Nz_)
                );
            } else {
                z_grid_ = spectral::SpectralDiffer::build_legendre_gauss_lobatto_nodes(
                        static_cast<int>(Nz_)
                );
            }
        }

        [[nodiscard]] std::vector<StateVec>
        solveSingleZ(size_t iz, const StateVec& y0) const {
            if (y0.empty()) {
                throw std::runtime_error("ZSliceSolver::solveSingleZ: empty initial state.");
            }

            ode::LHSFunc lhs = [this, iz](const StateVec& y, Real r, Real /*z_unused*/) {
                return (*builder_ptr_)(y, r, z_grid_[iz]);
            };

            const size_t nfields = y0.size();

            ode::SourceFunc src = [this, iz, nfields](Real r, Real z) {
                if (!source_func_) {
                    return StateVec(nfields, 0.0);
                }
                return source_func_(r, z_grid_[iz]);
            };

            ode::TransportODESystem sys{lhs, src};

            return ode::solve_single_z(r_grid_,
                                       y0,
                                       sys,
                                       z_grid_[iz],
                                       abs_tol_,
                                       rel_tol_,
                                       initial_step_);
        }

        const KerrMetricOutgoing& metric_;
        const RVector& r_grid_;
        size_t Nz_;
        const ghp::HeldBackgroundFieldsVectorized<CoordType>& held_;
        bool use_chebyshev_{false};

        RVector z_grid_{};

        const BuilderBase* builder_ptr_{nullptr};
        ode::SourceFunc source_func_{nullptr};

        Real abs_tol_{1e-12};
        Real rel_tol_{1e-12};
        Real initial_step_{1e-4};
    };

} // namespace ghz

#endif // HS1DATA_DATA_ZSLICESOLVER_HPP