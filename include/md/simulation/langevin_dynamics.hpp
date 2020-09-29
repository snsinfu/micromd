// Copyright snsinfu 2020.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SIMULATION_LANGEVIN_DYMNAMICS_HPP
#define MD_SIMULATION_LANGEVIN_DYMNAMICS_HPP

// This module provides a function for simulating Langevin dynamics.

#include <cstdint>
#include <functional>
#include <vector>

#include "../basic_types.hpp"
#include "../system.hpp"


namespace md
{
    // Parameters for Langevin dynamics simulation.
    struct langevin_dynamics_config
    {
        // Temperature of the environment in the unit of (thermal) energy. It
        // can be zero.
        md::scalar temperature = 1;

        // Time discretization step.
        md::scalar timestep = 1;

        // Number of steps to simulate.
        md::step steps = 1;

        // Seed for pseudo-random number generator.
        std::uint64_t seed = 0;

        // Optional callback function called after each step.
        std::function<void(md::step)> callback = {};
    };

    // simulate_langevin_dynamics simulates Langevin dynamics of the system.
    // It uses mass, friction, position and velocity particle attributes.
    inline void simulate_langevin_dynamics(
        md::system& system,
        md::langevin_dynamics_config config
    )
    {
        md::array_view<md::scalar const> masses = system.view_masses();
        md::array_view<md::scalar const> frictions = system.view_frictions();
        md::array_view<md::point> positions = system.view_positions();
        md::array_view<md::vector> velocities = system.view_velocities();

        md::scalar const temperature = config.temperature;
        md::scalar const timestep = config.timestep;
        md::index const particles = system.particle_count();

        std::vector<md::vector> forces(particles);
        md::random_engine random(config.seed);
        md::normal_distribution<md::scalar> normal;

        // BAOAB scheme.

        system.compute_force(forces);

        for (md::step step_ctr = 1; step_ctr <= config.steps; step_ctr++) {
            for (md::index i = 0; i < particles; i++) {
                md::scalar const damping = std::exp(-frictions[i] * timestep);
                md::scalar const agitation = 1 - damping * damping;
                md::scalar const sigma = std::sqrt(temperature * agitation / masses[i]);
                md::vector const normal_vector = {
                    normal(random), normal(random), normal(random)
                };

                // B step
                velocities[i] += 0.5 * timestep / masses[i] * forces[i];

                // A step
                positions[i] += 0.5 * timestep * velocities[i];

                // O step
                velocities[i] *= damping;
                velocities[i] += sigma * normal_vector;

                // A step
                positions[i] += 0.5 * timestep * velocities[i];
            }

            system.compute_force(forces);

            for (md::index i = 0; i < particles; i++) {
                // B step
                velocities[i] += 0.5 * timestep / masses[i] * forces[i];
            }

            if (config.callback) {
                config.callback(step_ctr);
            }
        }
    }
}

#endif
