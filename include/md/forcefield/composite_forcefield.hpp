// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_COMPOSITE_FORCEFIELD_HPP
#define MD_FORCEFIELD_COMPOSITE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that statically
// combines multiple forcefields into a single sum forcefield.

#include <iterator>
#include <numeric>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"


namespace md
{
    // composite_forcefield implements md::forcefield as a sum of zero or more
    // forcefield implementations given by the template parameter.
    //
    // Due to diamond inheritance the components must derive md::forcefield with
    // `virtual` keyword.
    template<typename... Components>
    class composite_forcefield : public virtual md::forcefield, public Components...
    {
    public:
        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            (void) system;

            md::scalar energies[] = {
                0,
                Components::compute_energy(system)...
            };
            return std::accumulate(std::begin(energies), std::end(energies), md::scalar(0));
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            int dummy[] = {
                0,
                (Components::compute_force(system, forces), 0)...
            };
            (void) dummy;
            (void) system;
            (void) forces;
        }
    };
}

#endif
