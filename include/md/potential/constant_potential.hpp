// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_CONSTANT_POTENTIAL_HPP
#define MD_POTENTIAL_CONSTANT_POTENTIAL_HPP

#include "../basic_types.hpp"


namespace md
{
    struct constant_potential
    {
        md::scalar energy = 0;

        md::scalar evaluate_energy(md::vector) const
        {
            return energy;
        }

        md::vector evaluate_force(md::vector) const
        {
            return {};
        }
    };
}

#endif