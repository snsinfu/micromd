// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_COMPOSITE_FORCEFIELD_HPP
#define MD_FORCEFIELD_COMPOSITE_FORCEFIELD_HPP

#include <iterator>
#include <numeric>

#include "../forcefield.hpp"
#include "../system.hpp"
#include "../typedef.hpp"


namespace md
{
    // composite_forcefield implements md::forcefield as a sum of zero or more
    // forcefield implementations given by the template parameter.
    //
    // Due to diamond inheritance the components must derive md::forcefield with
    // `virtual` keyword.
    template<typename... Components>
    class composite_forcefield : public virtual md::forcefield, Components...
    {
    public:
        // compute_energy computes the sum of the energy values computed by
        // child components.
        md::scalar compute_energy(md::system const& system) override
        {
            md::scalar energies[] = {
                0,
                Components::compute_energy(system)...
            };
            return std::accumulate(std::begin(energies), std::end(energies), md::scalar(0));
        }

        // compute_force computes the sum of the force values computed by child
        // Components.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            int dummy[] = {
                0,
                (Components::compute_force(system, forces), 0)...
            };
            (void) dummy;
        }
    };
}

#endif