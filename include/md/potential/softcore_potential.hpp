// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_SOFTCORE_POTENTIAL_HPP
#define MD_POTENTIAL_SOFTCORE_POTENTIAL_HPP

// This module provides a softcore short-range potential.

#include "../basic_types.hpp"
#include "../misc/math.hpp"


namespace md
{
    // softcore_potential implements the following bell-shaped short-range
    // potential energy function:
    //
    //     u(r) = e ( 1 - (r/s)^2 )^N               (r < s)
    //     F(r) = 2Ne/s^2 ( 1 - (r/s)^2 )^(N-1) r   (r < s)
    //
    // The potential energy is constantly zero at and beyond r = s.
    //
    // This function approximates the gaussian function
    //
    //    u(r) = e exp(- r^2 / (2 sigma^2) )
    //
    // when s = sqrt(2N) sigma and N is large.
    template<int N>
    struct softcore_potential
    {
        static_assert(N >= 1, "softcore_potential exponent must be positive");

        // The energy parameter e.
        md::scalar overlap_energy = 1;

        // The distance parameter s.
        md::scalar cutoff_distance = 1;


        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const k2 = 1 / (cutoff_distance * cutoff_distance);
            md::scalar const u2 = k2 * r.squared_norm();
            md::scalar const g = 1 - u2;

            if (g < 0) {
                return 0;
            }

            return overlap_energy * md::power<N>(g);
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const k2 = 1 / (cutoff_distance * cutoff_distance);
            md::scalar const u2 = k2 * r.squared_norm();
            md::scalar const g = 1 - u2;

            if (g < 0) {
                return {};
            }

            return 2 * N * overlap_energy * k2 * md::power<N - 1>(g) * r;
        }
    };
}

#endif
