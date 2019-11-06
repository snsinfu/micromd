// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_WCA_POTENTIAL_HPP
#define MD_POTENTIAL_WCA_POTENTIAL_HPP

// This module provides the Lennard-Jones potential.

#include "../basic_types.hpp"


namespace md
{
    // wca_potential implements the (12-6) Lennard-Jones potential function
    // cut off at the sigma:
    //
    //     u(r) = e ( (s/r)^12 - 2 (s/r)^6 + 1 )        (r < s),
    //     F(r) = 12 e ( (s/r)^12 - (s/r)^6 ) r / r^2   (r < s).
    //
    // The potential energy takes its minimum value 0 at r = s and beyond. This
    // potential energy is purely repulsive.
    struct wca_potential
    {
        // The enenrgy parameter e.
        md::scalar epsilon = 1;

        // The distance parameter s.
        md::scalar sigma = 1;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const u2 = sigma * sigma / r.squared_norm();
            md::scalar const u6 = u2 * u2 * u2;
            md::scalar const mod = u6 - 1;
            if (mod < 0) {
                return 0;
            }
            return epsilon * mod * mod;
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const r2_inv = 1 / r.squared_norm();
            md::scalar const u2 = sigma * sigma * r2_inv;
            md::scalar const u6 = u2 * u2 * u2;
            md::scalar const mod = u6 * u6 - u6;
            if (mod < 0) {
                return {};
            }
            return 12 * epsilon * mod * r2_inv * r;
        }
    };
}

#endif
