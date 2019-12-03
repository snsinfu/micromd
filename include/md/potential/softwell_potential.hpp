// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_SOFTWELL_POTENTIAL_HPP
#define MD_POTENTIAL_SOFTWELL_POTENTIAL_HPP

// This module provides a long-range well potential.

#include "../basic_types.hpp"
#include "../misc/math.hpp"


namespace md
{
    // `softwell_potential` is a bounded, long-range, attractive pairwise
    // potential:
    //
    //     u(r) = -e / (1 + (r/s)^p) .
    //
    // Near the zero distance this potential asymptotically behaves like the
    // p-th order harmonic potential. The potential decays to zero as the
    // distance `r` gets large.
    //
    template<int P = 2>
    struct softwell_potential
    {
        static_assert(P >= 2, "P must be >= 2");

        md::scalar energy = 1;
        md::scalar decay_distance = 1;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const k2 = 1 / (decay_distance * decay_distance);
            md::scalar const u2 = k2 * r.squared_norm();
            md::scalar const up = md::power_sqrt<P>(u2);
            md::scalar const v1 = 1 + up;
            return -energy / v1;
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const k2 = 1 / (decay_distance * decay_distance);
            md::scalar const u2 = k2 * r.squared_norm();
            md::scalar const up_2 = md::power_sqrt<P - 2>(u2);
            md::scalar const up = up_2 * u2;
            md::scalar const v1 = 1 + up;
            return (-P * energy * k2) * up_2 / (v1 * v1) * r;
        }
    };
}

#endif
