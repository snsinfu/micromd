// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_SEMISPRING_POTENTIAL_HPP
#define MD_POTENTIAL_SEMISPRING_POTENTIAL_HPP

// This module provides harmonic potential with near-side cutoff at the
// equilibrium distance.

#include "../basic_types.hpp"


namespace md
{
    // semispring_potential implements harmonic spring potential function with
    // near-side cutoff:
    //
    //     u(r) = K/2 (r - b)^2     (r > b),
    //     F(r) = -K (r - b)/r r    (r > b).
    //
    // b is the equilibrium distance.
    struct semispring_potential
    {
        md::scalar spring_constant = 1;
        md::scalar equilibrium_distance = 0;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const u = r.norm() - equilibrium_distance;
            if (u <= 0) {
                return 0;
            }
            return 0.5 * spring_constant * u * u;
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const r1 = r.norm();
            if (r1 <= equilibrium_distance) {
                return md::vector{};
            }
            return (spring_constant * equilibrium_distance / r1 - spring_constant) * r;
        }
    };
}

#endif
