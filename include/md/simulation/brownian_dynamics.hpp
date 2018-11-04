// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SIMULATION_BROWNIAN_DYMNAMICS_HPP
#define MD_SIMULATION_BROWNIAN_DYMNAMICS_HPP

#include <cstdint>
#include <functional>
#include <random>

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
}

#endif
