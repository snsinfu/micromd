// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_DIFF_POTENTIAL_HPP
#define MD_POTENTIAL_DIFF_POTENTIAL_HPP

// This module provides potential wrapper that computes difference of two
// potentials. Including this header enables binary `operator-` on potential
// functors.

#include "../basic_types.hpp"
#include "detail/detection.hpp"


namespace md
{
    // `diff_potential` is a potential functor that computes the difference of
    // two potential functors.
    template<typename Pot1, typename Pot2>
    struct diff_potential
    {
        Pot1 potential_1;
        Pot2 potential_2;

        diff_potential() = default;

        diff_potential(Pot1 const& pot1, Pot2 const& pot2)
            : potential_1{pot1}, potential_2{pot2}
        {
        }

        md::scalar evaluate_energy(md::vector r) const
        {
            return potential_1.evaluate_energy(r) - potential_2.evaluate_energy(r);
        }

        md::vector evaluate_force(md::vector r) const
        {
            return potential_1.evaluate_force(r) - potential_2.evaluate_force(r);
        }
    };

    // Make pairwise potential functors subtractable.
    template<typename Pot1, typename Pot2>
    diff_potential<
        detail::pairwise_potential_t<Pot1>,
        detail::pairwise_potential_t<Pot2>
    >
    operator-(Pot1 const& pot1, Pot2 const& pot2)
    {
        return diff_potential<Pot1, Pot2>{pot1, pot2};
    }
}

#endif
