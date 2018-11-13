// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SIMULATION_DETAIL_BROWNIAN_TIMESTEPPER_HPP
#define MD_SIMULATION_DETAIL_BROWNIAN_TIMESTEPPER_HPP

// This module provides timestepping algorithms used to implement for Brownian
// dynamics simulation.

#include <algorithm>
#include <cmath>

#include "../../basic_types.hpp"


namespace md
{
    namespace detail
    {
        // brownian_timestepper abstracts the algorithm for determining the
        // timestep of Brownian dynamics simulation.
        class brownian_timestepper
        {
        public:
            virtual ~brownian_timestepper() = default;

            // determine_timestep computes the timestep for given context.
            //
            // Params:
            //
            //   mobilities  = Array of mobilities of particles
            //   forces      = Array of forces acting on particles
            //   temperature = Temperature in the unit of energy
            //
            virtual md::scalar determine_timestep(
                md::array_view<md::scalar const> mobilities,
                md::array_view<md::vector const> forces,
                md::scalar temperature
            ) const = 0;
        };


        // monotonic_brownian_timestepper is a brownian_timestepper with a fixed
        // timestep.
        class monotonic_brownian_timestepper : public detail::brownian_timestepper
        {
        public:
            // Constructor takes the fixed timestep value.
            explicit monotonic_brownian_timestepper(md::scalar timestep)
                : timestep_{timestep}
            {
            }

            // determine_timestep returns the preconfigured timestep value.
            md::scalar determine_timestep(
                md::array_view<md::scalar const>,
                md::array_view<md::vector const>,
                md::scalar
            ) const override
            {
                return timestep_;
            }

        private:
            md::scalar timestep_;
        };


        // solve_brownian_timestep estimates the amount of time for a Brownian
        // particle to travel given distance. This function is used to implement
        // adaptive_brownian_timestepper below.
        //
        // Params:
        //
        //   distance    = Distance hte particle would travel
        //   mobility    = Mobility of the particle
        //   force2      = Squared force acting on the particle
        //   temperature = Temperature in the unit of energy
        //
        // Returns:
        //
        //   The estimated amount of time.
        //
        inline md::scalar solve_brownian_timestep(
            md::scalar distance,
            md::scalar mobility,
            md::scalar force2,
            md::scalar temperature
        )
        {
            constexpr md::scalar epsilon = 1e-6;
            constexpr md::scalar random_walk_factor = 2.55;

            md::scalar const a = mobility * mobility * force2;
            md::scalar const b = random_walk_factor * mobility * temperature;
            md::scalar const c = distance * distance;

            if (a * c < epsilon * b * b) {
                return c / b;
            }

            return (-b + std::sqrt(b * b + 4 * a * c)) / (2 * a);
        }


        // adaptive_brownian_timestepper changes timestep to limit the expected
        // particle dispalcement to specified maximum.
        class adaptive_brownian_timestepper : public detail::brownian_timestepper
        {
        public:
            // Constructor takes the upper bound of timestep and displacement.
            explicit adaptive_brownian_timestepper(md::scalar max_timestep, md::scalar distance)
                : max_timestep_{max_timestep}
                , distance_{distance}
            {
            }

            // determine_timestep returns timestep value that is small enough to
            // limit the expected displacement of particles within the upper
            // bound.
            md::scalar determine_timestep(
                md::array_view<md::scalar const> mobilities,
                md::array_view<md::vector const> forces,
                md::scalar temperature
            ) const override
            {
                md::scalar timestep = max_timestep_;

                for (md::index i = 0; i < mobilities.size(); i++) {
                    md::scalar dt = max_timestep_;

                    md::scalar const mobility = mobilities[i];
                    md::scalar const force2 = forces[i].squared_norm();

                    if (mobility != 0 || force2 != 0) {
                        dt = detail::solve_brownian_timestep(
                            distance_,
                            mobility,
                            force2,
                            temperature
                        );
                    }

                    timestep = std::min(timestep, dt);
                }

                return timestep;
            }

        private:
            md::scalar max_timestep_;
            md::scalar distance_;
        };
    }
}

#endif
