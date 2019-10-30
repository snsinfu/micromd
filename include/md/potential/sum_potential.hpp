// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_SUM_POTENTIAL_HPP
#define MD_POTENTIAL_SUM_POTENTIAL_HPP

// This module provides potential wrapper that computes sum of two potentials.
// Including this header enables `operator+` on potential functors.

#include "../basic_types.hpp"


namespace md
{
    // `sum_potential` is a potential functor that computes the sum of two
    // potential functors.
    template<typename Pot1, typename Pot2>
    struct sum_potential
    {
        Pot1 potential_1;
        Pot2 potential_2;

        sum_potential() = default;

        sum_potential(Pot1 const& pot1, Pot2 const& pot2)
            : potential_1{pot1}, potential_2{pot2}
        {
        }

        md::scalar evaluate_energy(md::vector r) const
        {
            return potential_1.evaluate_energy(r) + potential_2.evaluate_energy(r);
        }

        md::vector evaluate_force(md::vector r) const
        {
            return potential_1.evaluate_force(r) + potential_2.evaluate_force(r);
        }
    };

    namespace detail
    {
        // SFINAE utility that detects pairwise potential functor.
        template<
            typename Pot,
            md::scalar(Pot::*)(md::vector) const = &Pot::evaluate_energy,
            md::vector(Pot::*)(md::vector) const = &Pot::evaluate_force
        >
        Pot detect_pairwise_potential();

        template<typename Pot>
        using pairwise_potential_t = decltype(detect_pairwise_potential<Pot>());
    }

    // Make pairwise potential functors addable.
    template<typename Pot1, typename Pot2>
    sum_potential<
        detail::pairwise_potential_t<Pot1>,
        detail::pairwise_potential_t<Pot2>
    >
    operator+(Pot1 const& pot1, Pot2 const& pot2)
    {
        return sum_potential<Pot1, Pot2>{pot1, pot2};
    }
}

#endif
