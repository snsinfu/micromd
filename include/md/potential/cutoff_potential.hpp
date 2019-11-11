// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_CUTOFF_POTENTIAL_HPP
#define MD_POTENTIAL_CUTOFF_POTENTIAL_HPP

// This module provides potential wrapper that hard-cuts energy and force at a
// distance.

#include "../basic_types.hpp"


namespace md
{
    // `cutoff_potential` is a potential wrapper that cuts energy and force to
    // zero at `cutoff_distance`.
    template<typename Pot>
    struct cutoff_potential
    {
        Pot base;
        md::scalar cutoff_distance = 1;

        inline md::scalar evaluate_energy(md::vector r) const
        {
            if (should_cut(r)) {
                return 0;
            }
            return base.evaluate_energy(r);
        }

        inline md::vector evaluate_force(md::vector r) const
        {
            if (should_cut(r)) {
                return {};
            }
            return base.evaluate_force(r);
        }

    private:
        bool should_cut(md::vector r) const
        {
            return r.squared_norm() >= cutoff_distance * cutoff_distance;
        }
    };

    // `apply_cutoff` creates a `cutoff_potential` wrapping the given potential
    // object `pot`.
    template<typename Pot>
    md::cutoff_potential<Pot> apply_cutoff(Pot const& pot, md::scalar dcut)
    {
        return md::cutoff_potential<Pot>{pot, dcut};
    }
}

#endif
