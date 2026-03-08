//
// Created by Peter Zimmerman on 05.03.26.
//
#ifndef HS1DATA_DATA_BUILDERS_HPP
#define HS1DATA_DATA_BUILDERS_HPP
#pragma once


#include "ghz/transport/Operators.hpp"  // operators for the transport equation
#include "ghz/transport/Corrector.hpp"


#include <utility>
#include <vector>
#include <string>

namespace ghz {

/**
 *  @name findNearestRIndex
  * @brief Finds the nearest radial grid index for a given radial coordinate.
  * @param r The radial coordinate.
  * @return The index of the nearest radial grid point.
  */

    size_t findNearestRIndex(Real r, const RVector &r_grid) {
        if (r_grid.empty()) {
            throw std::invalid_argument("findNearestRIndex: r_grid is empty");
            // or: assert(!r_grid.empty());
        }
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
    size_t findNearestRIndexSorted(Real r, const RVector& r_grid) {
        if (r_grid.empty()) throw std::invalid_argument("r_grid empty");

        auto it = std::lower_bound(r_grid.begin(), r_grid.end(), r);
        if (it == r_grid.begin()) return 0;
        if (it == r_grid.end())   return r_grid.size() - 1;

        size_t hi = static_cast<size_t>(it - r_grid.begin());
        size_t lo = hi - 1;

        return (std::abs(r - r_grid[lo]) <= std::abs(r - r_grid[hi])) ? lo : hi; // tie -> lo
    }

/**
  * @class BuilderBase
  * @brief Abstract base class for building LHS operators for different components of the corrector.
  */

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

            dy[1] = (rho + rhob).real() * y[1] +
                    ((rho - rhob) * (rho - rhob)).real() * y[0]; // d^2 ReX/dr^2 = L[X_mmbar] applied to ReX
            dy[3] = (rho + rhob).real() * y[3] +
                    ((rho - rhob) * (rho - rhob)).real() * y[2]; // d^2 ImX/dr^2 = L[X_mmbar] applied to ImX

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
        BuilderXnm(const GHPSpectral &NX_mmbar,
                   const RVector &r_grid,
                   const RVector &z_grid)
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
  * This class computes lower order terms and source for the X_nn eqn.  \n
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
        const RVector r_grid_, z_grid_;

    }; // class BuilderXnn

}
#endif //HS1DATA_DATA_BUILDERS_HPP
