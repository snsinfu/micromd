// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_HEURISTICS_HPP
#define MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_HEURISTICS_HPP

// This module provides some heuristic helper functions for neighbor_list and
// subsystem_neighbor_list implementation.

#include "../../basic_types.hpp"
#include "../../misc/linear_hash.hpp"


namespace md
{
    namespace detail
    {
        // determine_hash returns a linear_hash object that is heuristically
        // parameterized to make neighbor_searcher perform good on given points.
        inline md::linear_hash determine_hash(
            md::array_view<md::point const> points,
            md::scalar dcut
        )
        {
            md::linear_hash hash;

            (void) dcut;

            // Benchmark simulations run fastest with this simple heuristic
            // among those I have tried.
            hash.modulus = md::linear_hash::uint(points.size() * 2 / 11);
            hash.modulus |= 1;

            return hash;
        }

        // determine_verlet_radius determines a good vertlet cutoff radius for
        // neighbor list calculation.
        inline md::scalar determine_verlet_radius(md::scalar dcut)
        {
            // Heuristic: 1.2 gives fairly good performance.
            return 1.2 * dcut;
        }
    }
}

#endif
