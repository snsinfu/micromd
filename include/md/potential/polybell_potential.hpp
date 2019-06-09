// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_POLYBELL_POTENTIAL_HPP
#define MD_POTENTIAL_POLYBELL_POTENTIAL_HPP

// This module provides a bell-shaped polynomial potential.

#include "../basic_types.hpp"
#include "../misc/math.hpp"


namespace md
{
    // Bell-shaped short-range potential:
    //
    //     u(r) = e ( 1 - (r/s)^P )^Q       (r < s)
    //
    // The potential energy is constantly zero at and beyond r = s. Larger P
    // makes the bell curve fatter, and larger Q makes the cutoff smoother.
    template<int P, int Q>
    struct polybell_potential
    {
        static_assert(P >= 2, "polybell exponent P must be >= 2");
        static_assert(Q >= 1, "polybell exponent Q must be positive");

        // The energy parameter e.
        md::scalar overlap_energy = 1;

        // The distance parameter s.
        md::scalar cutoff_distance = 1;


        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const k2 = 1 / (cutoff_distance * cutoff_distance);
            md::scalar const u2 = k2 * r.squared_norm();
            md::scalar const v = md::power_sqrt<P>(u2);
            md::scalar const g = 1 - v;

            if (g < 0) {
                return 0;
            }

            return overlap_energy * md::power<Q>(g);
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const k2 = 1 / (cutoff_distance * cutoff_distance);
            md::scalar const u2 = k2 * r.squared_norm();
            md::scalar const v = md::power_sqrt<P - 2>(u2);
            md::scalar const g = 1 - v * u2;

            if (g < 0) {
                return {};
            }

            return P * Q * overlap_energy * k2 * md::power<Q - 1>(g) * v * r;
        }
    };
}

#endif
