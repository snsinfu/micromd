// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_POWER_LAW_POTENTIAL_HPP
#define MD_POTENTIAL_POWER_LAW_POTENTIAL_HPP

#include "../basic_types.hpp"


namespace md
{
    namespace detail
    {
        inline md::scalar ramp(md::scalar x)
        {
            return x > 0 ? x : 0;
        }

        inline md::scalar power(md::scalar x, int n)
        {
            if (n == 0) {
                return 1;
            }

            md::scalar pow = x;
            for (int i = 1; i < n; i++) {
                pow *= x;
            }
            return pow;
        }
    }

    template<int N>
    struct power_law_potential
    {
        static_assert(N >= 1, "power_law_potential exponent must be positive");

        md::scalar overlap_energy = 1;
        md::scalar cutoff_distance = 1;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const k2 = 1 / (cutoff_distance * cutoff_distance);
            md::scalar const u2 = k2 * r.squared_norm();
            md::scalar const g = 1 - u2;

            if (g < 0) {
                return 0;
            }

            return overlap_energy * detail::power(g, N);
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const k2 = 1 / (cutoff_distance * cutoff_distance);
            md::scalar const u2 = k2 * r.squared_norm();
            md::scalar const g = 1 - u2;

            if (g < 0) {
                return {};
            }

            return N * overlap_energy * k2 * detail::power(g, N - 1) * r;
        }
    };
}

#endif
