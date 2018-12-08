// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_SPRING_POTENTIAL_HPP
#define MD_POTENTIAL_SPRING_POTENTIAL_HPP

// This module provides harmonic potential with nonzero equilibriums distance.

#include "../basic_types.hpp"


namespace md
{
    // spring_potential implements harmonic spring potential function:
    //
    //     u(r) = K/2 (r - b)^2 ,
    //     F(r) = -K (r - b)/r r .
    //
    // b is the equilibrium distance.
    struct spring_potential
    {
        md::scalar spring_constant = 1;
        md::scalar equilibrium_distance = 0;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const u = r.norm() - equilibrium_distance;
            return 0.5 * spring_constant * u * u;
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const r1 = r.norm();
            if (r1 == 0) {
                return md::vector{};
            }
            return (spring_constant * equilibrium_distance / r1 - spring_constant) * r;
        }
    };
}

#endif
