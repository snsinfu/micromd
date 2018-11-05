// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SIMULATION_BROWNIAN_DYMNAMICS_HPP
#define MD_SIMULATION_BROWNIAN_DYMNAMICS_HPP

#include <cmath>
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
        // Temperature of the environment in the unit of (thermal) energy. It
        // can be zero.
        md::scalar temperature = 1;

        // Time discretization step.
        md::scalar timestep = 1;

        // Set this to nonzero to adaptively determine timestep so that the mean
        // displacement is this spacestep value.
        md::scalar spacestep = 0;

        // Number of steps to simulate.
        md::step steps = 1;

        // Seed for pseudo-random number generator.
        std::uint64_t seed = 0;

        // Optional function called after each step.
        std::function<void(md::step)> callback;
    };

    namespace detail
    {
        inline md::scalar solve_brownian_timestep(
            md::scalar spacestep,
            md::scalar mobility,
            md::scalar force2,
            md::scalar temperature
        )
        {
            md::scalar const a = mobility * mobility * force2;
            md::scalar const b = 2.55 * mobility * temperature;
            md::scalar const c = spacestep * spacestep;

            md::scalar const ac_threshold = 1e-6 * b * b;
            if (a * c < ac_threshold) {
                return c / b;
            }

            return (-b + std::sqrt(b * b + 4 * a * c)) / (2 * a);
        }

        inline md::scalar determine_brownian_timestep(
            md::scalar max_timestep,
            md::scalar spacestep,
            md::array_view<md::scalar const> mobilities,
            md::array_view<md::vector const> forces,
            md::scalar temperature
        )
        {
            if (spacestep == 0) {
                return max_timestep;
            }

            md::scalar timestep = max_timestep;

            for (md::index i = 0; i < mobilities.size(); i++) {
                md::scalar dt = max_timestep;

                md::scalar const mobility = mobilities[i];
                md::scalar const force2 = forces[i].squared_norm();

                if (mobility != 0 || force2 != 0) {
                    dt = detail::solve_brownian_timestep(spacestep, mobility, force2, temperature);
                }

                if (dt < timestep) {
                    timestep = dt;
                }
            }

            return timestep;
        }
    }

    class brownian_dynamics_simulator
    {
    public:
        explicit brownian_dynamics_simulator(md::brownian_dynamics_config config)
            : config_{config}
        {
        }

        void run_simulation(md::system& system)
        {
            forces_ = std::vector<md::vector>(system.particle_count());
            weiners_ = std::vector<md::vector>(system.particle_count());

            for (md::step step = 0; step < config_.steps; step++) {
                simulate_step(system);
                callback_step(step);
            }
        }

    private:
        void simulate_step(md::system& system)
        {
            md::array_view<md::scalar const> mobilities = system.view_mobilities();
            md::array_view<md::point> positions = system.view_positions();

            system.compute_force(forces_);

            md::scalar const timestep = detail::determine_brownian_timestep(
                config_.timestep,
                config_.spacestep,
                mobilities,
                forces_,
                config_.temperature
            );

            std::normal_distribution<md::scalar> normal;

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::scalar const mu_dt = timestep * mobilities[i];
                md::scalar const sigma = std::sqrt(2 * config_.temperature * mu_dt);
                md::vector const weiner = {
                    sigma * normal(random_),
                    sigma * normal(random_),
                    sigma * normal(random_)
                };

                positions[i] += mu_dt * forces_[i];
                positions[i] += 0.5 * (weiner + weiners_[i]);
                weiners_[i] = weiner;
            }
        }

        void callback_step(md::step step) const
        {
            if (config_.callback) {
                config_.callback(step);
            }
        }

    private:
        md::brownian_dynamics_config config_;
        std::mt19937_64 random_;
        std::vector<md::vector> forces_;
        std::vector<md::vector> weiners_;
    };

    // simulate_brownian_dynamics simulates Brownian dynamics of the system.
    inline void simulate_brownian_dynamics(md::system& system, md::brownian_dynamics_config config)
    {
        md::brownian_dynamics_simulator{config}.run_simulation(system);
    }
}

#endif
