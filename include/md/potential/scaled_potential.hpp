// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_SCALED_POTENTIAL_HPP
#define MD_POTENTIAL_SCALED_POTENTIAL_HPP

// This module provides potential wrapper that computes scalar-multiplied
// potential. Including this header enables `operator*` on potential functors.

#include "../basic_types.hpp"
#include "detail/detection.hpp"


namespace md
{
    // `scaled_potential` is a potential functor that computes scalar-multiplied
    // potential functor.
    template<typename Pot>
    struct scaled_potential
    {
        Pot potential;
        md::scalar factor;

        scaled_potential() = default;

        scaled_potential(Pot const& pot, md::scalar factor_)
            : potential{pot}, factor{factor_}
        {
        }

        md::scalar evaluate_energy(md::vector r) const
        {
            return factor * potential.evaluate_energy(r);
        }

        md::vector evaluate_force(md::vector r) const
        {
            return factor * potential.evaluate_force(r);
        }
    };

    // Make pairwise potential functors scalar-multipliable from right.
    template<typename Pot>
    scaled_potential<detail::pairwise_potential_t<Pot>>
    operator*(Pot const& pot, md::scalar factor)
    {
        return scaled_potential<Pot>{pot, factor};
    }

    // Make pairwise potential functors scalar-multipliable from left.
    template<typename Pot>
    scaled_potential<detail::pairwise_potential_t<Pot>>
    operator*(md::scalar factor, Pot const& pot)
    {
        return scaled_potential<Pot>{pot, factor};
    }

    // Make pairwise potential functors negatable.
    template<typename Pot>
    scaled_potential<detail::pairwise_potential_t<Pot>> operator-(Pot const& pot)
    {
        return scaled_potential<Pot>{pot, -1};
    }
}

#endif
