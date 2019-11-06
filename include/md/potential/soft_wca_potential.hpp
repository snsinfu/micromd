// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_SOFT_WCA_POTENTIAL_HPP
#define MD_POTENTIAL_SOFT_WCA_POTENTIAL_HPP

// This module provides a soft variation of the WCA potential.

#include "../basic_types.hpp"


namespace md
{
    // soft_wca_potential is a soft variation of the WCA potential function:
    //
    //     u(r) = e ((k + 1)/(k + (r/s)^6) - 1)^2 .
    //
    // The parameter k controls the soft-ness of the repulsion:
    //
    //     u(0) = e/k^2 .
    //
    // The potential matches the WCA potential when k=0. The potential is purely
    // repulsive and vanishes at and beyond r=s.
    struct soft_wca_potential
    {
        md::scalar epsilon = 1;
        md::scalar sigma = 1;
        md::scalar softness = 0.1;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const u2 = r.squared_norm() / (sigma * sigma);
            md::scalar const u6 = u2 * u2 * u2;
            md::scalar const g = (softness + 1) / (softness + u6) - 1;
            if (g < 0) {
                return 0;
            }
            return epsilon * g * g;
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const sigma_inv2 = 1 / (sigma * sigma);
            md::scalar const u2 = r.squared_norm() * sigma_inv2;
            md::scalar const u4 = u2 * u2;
            md::scalar const u6 = u4 * u2;
            md::scalar const g = (softness + 1) / (softness + u6);
            if (g * g < g) {
                return {};
            }
            return 12 * epsilon * sigma_inv2 * u4 * (g * g - g) * r;
        }
    };
}

#endif
