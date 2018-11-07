// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SIMULATION_DETAIL_BROWNIAN_SIMULATOR_HPP
#define MD_SIMULATION_DETAIL_BROWNIAN_SIMULATOR_HPP

// This module provides a helper class to simulate Brownian dynamics.
//
// detail::brownian_simulator

#include <cmath>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "../../basic_types.hpp"
#include "../../system.hpp"

#include "brownian_timestepper.hpp"


namespace md
{
    namespace detail
    {
        class brownian_simulator
        {
        public:
            // Params:
            //
            //   system      = Simulated system
            //   timestepper = Timestepper to use
            //   temperature = Temperature in the unit of energy
            //   seed        = Random seed
            //
            brownian_simulator(
                md::system& system,
                detail::brownian_timestepper& timestepper,
                md::scalar temperature,
                std::uint64_t seed
            )
                : system_{system}
                , timestepper_{timestepper}
                , temperature_{temperature}
            {
                std::seed_seq seed_seq {
                    std::uint32_t(seed >> 32),
                    std::uint32_t(seed)
                };
                random_.seed(seed_seq);

                forces_.resize(system.particle_count());
                weiners_.resize(system.particle_count());
            }

            // simulate_step simulates a discretized step.
            void simulate_step()
            {
                md::array_view<md::scalar const> mobilities = system_.view_mobilities();
                md::array_view<md::point> positions = system_.view_positions();
                md::array_view<md::vector> forces = forces_;
                md::array_view<md::vector> weiners = weiners_;

                // Second-order BAOAB limit scheme.

                system_.compute_force(forces);

                md::scalar const timestep = timestepper_.determine_timestep(
                    mobilities,
                    forces,
                    temperature_
                );

                std::normal_distribution<md::scalar> normal;

                for (md::index i = 0; i < system_.particle_count(); i++) {
                    md::scalar const mu_dt = timestep * mobilities[i];
                    md::scalar const sigma = std::sqrt(2 * temperature_ * mu_dt);
                    md::vector const weiner = {
                        sigma * normal(random_),
                        sigma * normal(random_),
                        sigma * normal(random_)
                    };

                    positions[i] += mu_dt * forces[i];
                    positions[i] += 0.5 * (weiner + weiners[i]);
                    weiners[i] = weiner;
                }
            }

        private:
            md::system& system_;
            detail::brownian_timestepper& timestepper_;
            md::scalar temperature_;
            std::vector<md::vector> forces_;
            std::vector<md::vector> weiners_;
            std::mt19937_64 random_;
        };
    }
}

#endif