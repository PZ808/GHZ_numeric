//
// Created by Peter Zimmerman on 05.03.26.
//

#ifndef HS1DATA_DATA_ZSLICESOLVER_HPP
#define HS1DATA_DATA_ZSLICESOLVER_HPP

#include "ghz/transport/Corrector.hpp"
#include "ghz/transport/Builders.hpp"

namespace ghz {
    using namespace teuk::literals;
    using StateVec = ode::StateVec;
    using StateMat = std::vector <std::vector<StateVec>>;
    using RVector = std::vector<Real>;
    using CVector = std::vector<Complex>;
    using GHPSpectral = spectral::SpectralGHPVectorized;

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
                     const ghp::HeldBackgroundFieldsVectorized <CoordType> &held,
                     bool use_chebyshev = false)
                : metric_(metric), r_grid_(r_grid), Nz_(Nz), held_(held), use_chebyshev_(use_chebyshev) {
            buildZGrid();
        }

        StateMat solve() {
            StateMat sol_2D(Nz_);

#pragma omp parallel for
            for (size_t iz = 0; iz < Nz_; ++iz) {
                sol_2D[iz] = solveSingleZ(iz);
            }
            return sol_2D;
        }

        void setInitialCondition(const ode::StateVec &y0) {
            y0_ = y0;
            Nfields_ = static_cast<int>(y0_.size());
        }

        void setSourceData(const std::vector <std::vector<Real>> &source) { source_data_ = source; }

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

        std::vector <ode::StateVec> solveSingleZ(size_t iz) {
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
        const std::vector <Real> &r_grid_;
        size_t Nz_;
        bool use_chebyshev_ = false;
        const ghp::HeldBackgroundFieldsVectorized <CoordType> &held_;
        std::vector <Real> z_grid_;
        const BuilderBase *builder_ptr_ = nullptr;
        ode::SourceFunc source_func_ = nullptr;
        std::vector <std::vector<Real>> source_data_;
        int Nfields_ = 4;
        ode::StateVec y0_;
        Real abs_tol_ = 1e-12;
        Real rel_tol_ = 1e-12;
        Real initial_step_ = 1e-4;
    }; // class ZSliceSolver
}
#endif //HS1DATA_DATA_ZSLICESOLVER_HPP
