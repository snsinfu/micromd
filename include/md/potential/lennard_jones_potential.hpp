// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_LENNARD_JONES_POTENTIAL_HPP
#define MD_POTENTIAL_LENNARD_JONES_POTENTIAL_HPP

// This module provides the Lennard-Jones potential.

#include "../basic_types.hpp"


namespace md
{
    // lennard_jones_potential implements the (12-6) Lennard-Jones potential
    // function:
    //
    //     u(r) = e ( (s/r)^12 - 2 (s/r)^6 ) ,
    //     F(r) = 12 e ( (s/r)^12 - (s/r)^6 ) r / r^2 .
    //
    // The potential energy takes its minimum value -e at r = s.
    struct lennard_jones_potential
    {
        // The enenrgy parameter e.
        md::scalar epsilon = 1;

        // The distance parameter s.
        md::scalar sigma = 1;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const u2 = sigma * sigma / r.squared_norm();
            md::scalar const u6 = u2 * u2 * u2;

            return epsilon * (u6 * u6 - u6 - u6);
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const r2_inv = 1 / r.squared_norm();
            md::scalar const u2 = sigma * sigma * r2_inv;
            md::scalar const u6 = u2 * u2 * u2;

            return 12 * epsilon * (u6 * u6 - u6) * r2_inv * r;
        }
    };
}

#endif
