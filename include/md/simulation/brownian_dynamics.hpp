// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SIMULATION_BROWNIAN_DYMNAMICS_HPP
#define MD_SIMULATION_BROWNIAN_DYMNAMICS_HPP

// This module provides a function for simulating Brownian dynamics.

#include <cstdint>
#include <functional>
#include <memory>

#include "../basic_types.hpp"
#include "../system.hpp"

#include "detail/brownian_simulator.hpp"
#include "detail/brownian_timestepper.hpp"


namespace md
{
    // brownian_dynamics_config holds Brownian dynamics parameters.
    struct brownian_dynamics_config
    {
        // Temperature of the environment in the unit of (thermal) energy. It
        // can be zero.
        md::scalar temperature = 1;

        // Time discretization step.
        md::scalar timestep = 1;

        // Set this to nonzero to adaptively determine timestep so that the mean
        // displacement is bounded by spacestep.
        md::scalar spacestep = 0;

        // Number of steps to simulate.
        md::step steps = 1;

        // Seed for pseudo-random number generator.
        std::uint64_t seed = 0;

        // Optional function called after each step.
        std::function<void(md::step)> callback = {};
    };

    // simulate_brownian_dynamics simulates Brownian dynamics of the system.
    inline void simulate_brownian_dynamics(md::system& system, md::brownian_dynamics_config config)
    {
        std::unique_ptr<detail::brownian_timestepper> timestepper;

        if (config.spacestep != 0) {
            timestepper = std::make_unique<detail::adaptive_brownian_timestepper>(
                config.timestep,
                config.spacestep
            );
        } else {
            timestepper = std::make_unique<detail::monotonic_brownian_timestepper>(
                config.timestep
            );
        }

        detail::brownian_simulator simulator(system, *timestepper, config.temperature, config.seed);

        for (md::step step_ctr = 0; step_ctr < config.steps; step_ctr++) {
            simulator.simulate_step();

            if (config.callback) {
                config.callback(step_ctr);
            }
        }
    }
}

#endif
