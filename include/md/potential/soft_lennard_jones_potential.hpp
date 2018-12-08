// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_SOFT_LENNARD_JONES_POTENTIAL_HPP
#define MD_POTENTIAL_SOFT_LENNARD_JONES_POTENTIAL_HPP

// This module provides a soft variation of the Lennard-Jones potential.

#include "../basic_types.hpp"


namespace md
{
    // soft_lennard_jones_potential is a soft variation of the (12-6) Lennard-
    // Jones potential function:
    //
    //     u(r) = e ((k + 1)/(k + (r/s)^6) - 1)^2 .
    //
    // The parameter k controls the soft-ness of the repulsion:
    //
    //     u(0) = e/k^2 .
    //
    // The potential matches the LJ potential when k=0.
    struct soft_lennard_jones_potential
    {
        md::scalar epsilon = 1;
        md::scalar sigma = 1;
        md::scalar softness = 0.1;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const u2 = r.squared_norm() / (sigma * sigma);
            md::scalar const u6 = u2 * u2 * u2;
            md::scalar const g = (softness + 1) / (softness + u6) - 1;
            return epsilon * g * g;
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const sigma_inv2 = 1 / (sigma * sigma);
            md::scalar const u2 = r.squared_norm() * sigma_inv2;
            md::scalar const u4 = u2 * u2;
            md::scalar const u6 = u4 * u2;
            md::scalar const g = (softness + 1) / (softness + u6);
            return 12 * epsilon * sigma_inv2 * u4 * (g * g - g) * r;
        }
    };
}

#endif
