// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SIMULATION_BROWNIAN_DYMNAMICS_HPP
#define MD_SIMULATION_BROWNIAN_DYMNAMICS_HPP

#include <cstdint>
#include <functional>
#include <random>
#include <vector>

#include "../basic_types.hpp"
#include "../system.hpp"


namespace md
{
    // brownian_dynamics_config holds Brownian dynamics parameters.
    struct brownian_dynamics_config
    {
        md::scalar temperature = 1;
        md::scalar timestep = 1;
        md::step steps = 1;
        std::uint64_t seed = 0;
        std::function<void(md::step)> callback;
    };

    // simulate_brownian_dynamics simulates Brownian dynamics of the system.
    inline void simulate_brownian_dynamics(md::system& system, md::brownian_dynamics_config config)
    {
        std::seed_seq seed {
            unsigned(config.seed >> 32),
            unsigned(config.seed)
        };
        std::mt19937_64 random(seed);
        std::normal_distribution<md::scalar> normal;

        // BAOAB limit scheme.

        md::array_view<md::point> positions = system.view_positions();
        md::array_view<md::scalar> mobilities = system.view_mobilities();
        std::vector<md::vector> forces(system.particle_count());
        std::vector<md::vector> weiners(system.particle_count());

        for (md::step step = 0; step < config.steps; step++) {
            system.compute_force(forces);

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::scalar const mu_dt = config.timestep * mobilities[i];
                md::scalar const sigma = std::sqrt(2 * config.temperature * mu_dt);
                md::vector const weiner = {
                    sigma * normal(random),
                    sigma * normal(random),
                    sigma * normal(random)
                };

                positions[i] += mu_dt * forces[i];
                positions[i] += 0.5 * (weiner + weiners[i]);
                weiners[i] = weiner;
            }

            if (config.callback) {
                config.callback(step);
            }
        }
    }
}

#endif
