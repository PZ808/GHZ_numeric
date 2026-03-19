//
// Created by Peter Zimmerman on 05.03.26.
//
#ifndef HS1DATA_DATA_BUILDERS_HPP
#define HS1DATA_DATA_BUILDERS_HPP
#pragma once


#include "Operators.hpp"  // operators for the transport equation
#include "ghz/transport/Corrector.hpp"


#include <utility>
#include <vector>
#include <string>
#include <algorithm>

namespace ghz {

/**
 *  @name findNearestRIndex
  * @brief Finds the nearest radial grid index for a given radial coordinate.
  * @param r The radial coordinate.
  * @return The index of the nearest radial grid point.
  */

    inline size_t findNearestRIndex(Real r, const RVector &r_grid) {
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
    inline size_t findNearestRIndexSorted(Real r, const RVector& r_grid) {
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
        explicit BuilderBase(Real a) : a_(a) {}
        virtual ~BuilderBase() = default;

        virtual StateVec operator()(const StateVec& y, Real r, Real z) const = 0;

    protected:
        // rho in kinnersley
        [[nodiscard]] Complex rho_K(Real r, Real z) const {
            return -Real(1) / (r - Complex(0.0, a_ * z));
        }

    private:
        Real a_;
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
        explicit BuilderXmmbar(Real a) : BuilderBase(a) {}
        StateVec operator()(const StateVec &y,
                            Real r, Real z) const override {

            StateVec dy(4, 0.0); // initially zero. fill with LHS deriv operator applied to y

            const Complex rho = rho_K(r, z);
            const Complex rhob = std::conj(rho);

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
   * This class computes the LHS operator for the X_nm level, for the eqn.  \n
   * rho/(2(rho + rhob)) * d/dr ( (rho + rhob)^2 * d/dr [X_nm/(rho*(rho + rhob))] ) = T_lm + N[x_mmbar] \n
   * L[f] = rho dX/dr + rho^2 f - rhob*(rho-rhob)f
   */
    class BuilderXnm : public BuilderBase {

    public:
        explicit BuilderXnm(Real a) : BuilderBase(a) {}

        ode::StateVec operator()(const StateVec &y,
                                 Real r, Real z) const override {
            // Expect y.size() == 4 (ReX, ReX', ImX, ImX')
            ode::StateVec dy(4, 0.0);

            // unpack X_nm into complex quantities
            Complex Xnm(y[0], y[2]); // ReXnm and ImXnm state vector
            Complex Vnm(y[1], y[3]); // dReX_nm/dr and dImXnm/dr
            const Complex rho = rho_K(r, z);
            const Complex rhob = std::conj(rho);


            Complex LXnm = 2.0_r * (
                    rho*rho * Xnm
                   - rhob*(rho-rhob) * Xnm
                   + rho * Vnm
                   );

            // pack dy
            dy[0] = std::real(Vnm);
            dy[2] = std::imag(Vnm);
            dy[1] = std::real(LXnm);
            dy[3] = std::imag(LXnm);

            return dy;
        }


    };// class BuilderXnm

/**
  * @class BuilderXnn
  * @brief Concrete implementation of BuilderBase for X_nn  \n
  *
  *  1/2 * (rho+rhob)^2 * Thorn[ X_nn/(rho+rhob) ] = Tln + Re U[X_mmbar] + Re V[X_nm] \n
  * L[f] = - rho*rhob/r^2 * ( r^2-a^2*z^2 ) * f
  */
    class BuilderXnn : public BuilderBase {
    public:
        explicit BuilderXnn(Real a) : BuilderBase(a) {}

        ode::StateVec operator()(const ode::StateVec &y,
                                 Real r, Real z) const override {

            // State: y = {Re Xnn, Im Xnn}
            Complex Xnn(y[0], y[1]);

            //Real Xnn = y[0];
            //ode::StateVec dy(1);

            const Complex rho = rho_K(r, z);
            const Complex rhob = std::conj(rho);

            // compute dXnn/dr from first-order formula
            Complex dX = ( (rho*rho + rhob*rhob)/(rho+rhob) ) * Xnn;
            //Real dX = std::real(((rho*rho + rhob*rhob)/(rho+rhob))) * Xnn;
            //dy[0] = dX;

            ode::StateVec dy(2);
            dy[0] = std::real(dX);
            dy[1] = std::imag(dX);
            return dy;
        }

    }; // class BuilderXnn

}
#endif //HS1DATA_DATA_BUILDERS_HPP
