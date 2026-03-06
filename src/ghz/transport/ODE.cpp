//
// Created by Peter Zimmerman on 06.03.26.
//

#include "ghz/transport/ODE.hpp"


namespace ode {

    // Adaptive solver for coupled system over r for a single z
    std::vector<StateVec> solve_single_z(
            const std::vector<Real>& r_grid, // nodes in r or (sigma) or the chebyshev grid relative to these values
            const StateVec& y0,              // initial conditions at r_min
            const TransportODESystem& ode,   // ODE system
            Real z,                          // fixed z value
            Real abs_tol=1e-9,             // solver abs tolerance
            Real rel_tol=1e-9,             // solver rel tolerance
            Real initial_step=1e-3) {      // initial step size

        using namespace boost::numeric::odeint;
        // Define the stepper type for adaptive integration
        using stepper_type = boost::numeric::odeint::runge_kutta_cash_karp54<
                StateVec, Real, StateVec, Real,
                boost::numeric::odeint::range_algebra,
                boost::numeric::odeint::default_operations>;

        StateVec y = y0;  // initialize state vector
        size_t Nr = r_grid.size(); // number of r grid points
        std::vector<StateVec> solution; // store solution at r grid points
        solution.reserve(Nr); // reserve space

        size_t idx = 0;
        // Observer lambda to capture solution at specified r_grid points
        auto observer = [&](const StateVec &y_val, Real r_val){
            while (idx < Nr && r_val >= r_grid[idx]) {
                solution.push_back(y_val);
                ++idx;
            }
        };
        // Wrap the ODE to fix z
        auto ode_wrapper = [&](const StateVec &y_in,
                               StateVec &dydr, Real r) {
            ode(y_in, dydr, r, z);
        };

        // Perform adaptive integration over r
        integrate_adaptive(make_controlled(abs_tol, rel_tol, stepper_type()),
                           ode_wrapper, y,
                           r_grid.front(), r_grid.back(),
                           initial_step,
                           observer);

        // wrap up: ensure solution has entries for all r_grid points
        while (solution.size() < Nr)
            solution.push_back(y); // fill with last value

        return solution;
    } // solve_single_z

    // Solver over r_grid and z_grid
    std::vector<std::vector<StateVec>> solve_2D(const std::vector<Real>& r_grid, // nodes in r
                                                const std::vector<Real>& z_grid, // nodes in z
                                                const StateVec& y0,
                                                const TransportODESystem& ode,
                                                Real abs_tol = 1e-9,
                                                Real rel_tol = 1e-9,
                                                Real initial_step = 1e-3) {
        size_t Nz = z_grid.size();
        // container for full 2D solution
        std::vector<std::vector<StateVec>> sol_2D; // sol_2D[iz][ir][component]
        sol_2D.reserve(Nz);

        // loop over z grid points
        for (size_t iz=0; iz<Nz; ++iz){
            Real z = z_grid[iz];
            sol_2D.push_back(solve_single_z(r_grid, y0, ode, z, abs_tol, rel_tol, initial_step));
        }
        return sol_2D; // sol_2D[iz][ir][component]
    }

}