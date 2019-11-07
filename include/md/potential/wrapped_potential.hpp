// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_WRAPPED_POTENTIAL_HPP
#define MD_POTENTIAL_WRAPPED_POTENTIAL_HPP

// This module provides identity potential wrapper.

#include "../basic_types.hpp"
#include "detail/detection.hpp"


namespace md
{
    // `wrapped_potential` is a potential functor that computes exactly the same
    // energy and force as the wrapped potential. Wrapping custom potential with
    // this class enables arithmetic operations just like other potentials in
    // the `md` namespace.
    template<typename Pot>
    struct wrapped_potential
    {
        Pot potential;

        wrapped_potential() = default;

        explicit wrapped_potential(Pot const& pot)
            : potential{pot}
        {
        }

        inline md::scalar evaluate_energy(md::vector r) const
        {
            return potential.evaluate_energy(r);
        }

        inline md::vector evaluate_force(md::vector r) const
        {
            return potential.evaluate_force(r);
        }
    };

    // Make `wrapped_potential`.
    template<typename Pot>
    wrapped_potential<Pot> wrap_potential(Pot const& pot)
    {
        return wrapped_potential<Pot>{pot};
    }
}

#endif
