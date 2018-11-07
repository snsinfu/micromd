// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SIMULATION_NEWTONIAN_DYMNAMICS_HPP
#define MD_SIMULATION_NEWTONIAN_DYMNAMICS_HPP

#include <cmath>
#include <functional>
#include <vector>

#include "../basic_types.hpp"
#include "../system.hpp"


namespace md
{
    // Parameters for Newtonian dynamics simulation.
    struct newtonian_dynamics_config
    {
        // Time discretization step.
        md::scalar timestep = 1;

        // Number of steps to simulate.
        md::step steps = 1;

        // Optional function called after each step.
        std::function<void(md::step)> callback;
    };

    inline void simulate_newtonian_dynamics(md::system& system, md::newtonian_dynamics_config config)
    {
        std::vector<md::vector> forces(system.particle_count());

        md::array_view<md::scalar const> masses = system.view_masses();
        md::array_view<md::point> positions = system.view_positions();
        md::array_view<md::vector> velocities = system.view_velocities();

        // Velocity Verlet scheme.

        for (md::step step_ctr = 0; step_ctr < config.steps; step_ctr++) {
            md::scalar const timestep = config.timestep;

            for (md::index i = 0; i < system.particle_count(); i++) {
                velocities[i] += timestep / (2 * masses[i]) * forces[i];
                positions[i] += timestep * velocities[i];
            }

            system.compute_force(forces);

            for (md::index i = 0; i < system.particle_count(); i++) {
                velocities[i] += timestep / (2 * masses[i]) * forces[i];
            }

            if (config.callback) {
                config.callback(step_ctr);
            }
        }
    }
}

#endif
