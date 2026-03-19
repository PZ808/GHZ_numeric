//
// Created by Peter Zimmerman on 06.03.26.
//

#include "ghz/transport/RK/ODE.hpp"


namespace ode {

    inline std::vector<StateVec> solve_single_iz(
            const std::vector<Real>& r_grid,
            const StateVec& y0,
            const TransportODESystem& ode,
            size_t iz,
            Real z,
            Real abs_tol,
            Real rel_tol,
            Real initial_step)
    {
        using namespace boost::numeric::odeint;

        if (r_grid.empty()) {
            throw std::invalid_argument("solve_single_z: r_grid is empty.");
        }
        if (y0.empty()) {
            throw std::invalid_argument("solve_single_z: y0 is empty.");
        }
        if (r_grid.size() == 1) {
            return {y0};
        }
        if (!std::is_sorted(r_grid.begin(), r_grid.end())) {
            throw std::invalid_argument("solve_single_z: r_grid must be sorted ascending.");
        }

        using stepper_type = runge_kutta_cash_karp54<
                StateVec, Real, StateVec, Real,
                range_algebra, default_operations>;

        auto controlled = make_controlled(
                abs_tol, rel_tol,
                stepper_type());

        StateVec y = y0;
        std::vector<StateVec> solution;
        solution.reserve(r_grid.size());
        solution.push_back(y0); // value at r_grid[0]

        auto ode_wrapper = [&](const StateVec& y_in, StateVec& dydr, Real r) {
            ode(y_in, dydr, r, z);
        };

        Real h = initial_step;

        for (size_t i=0; i<r_grid.size()-1; ++i) {
            const Real r0 = r_grid[i];
            const Real r1 = r_grid[i+1];

            integrate_adaptive(controlled, ode_wrapper, y, r0, r1, h);

            solution.push_back(y); // now this really is the state at r1
        }

        return solution;
    }// solve_single_z

}