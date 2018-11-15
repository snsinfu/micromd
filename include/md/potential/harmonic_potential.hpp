// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_HARMONIC_POTENTIAL_HPP
#define MD_POTENTIAL_HARMONIC_POTENTIAL_HPP

// This module provides the harmonic potential.

#include "../basic_types.hpp"


namespace md
{
    // harmonic_potential implements the centered harmonic potential function:
    //
    //     u(r) = K/2 r^2 ,
    //     F(r) = -K r .
    //
    struct harmonic_potential
    {
        // The spring constant K.
        md::scalar spring_constant = 1;

        md::scalar evaluate_energy(md::vector r) const
        {
            return 0.5 * spring_constant * r.squared_norm();
        }

        md::vector evaluate_force(md::vector r) const
        {
            return -spring_constant * r;
        }
    };
}

#endif
