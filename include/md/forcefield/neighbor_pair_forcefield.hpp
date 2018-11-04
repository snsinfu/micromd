// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_NEIGHBOR_PAIR_FORCEFIELD_HPP
#define MD_FORCEFIELD_NEIGHBOR_PAIR_FORCEFIELD_HPP

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/neighbor_list.hpp"


namespace md
{
    template<typename Derived>
    class neighbor_pair_forcefield : public virtual md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view(md::position_attribute);
            md::scalar sum = 0;

            for (auto const pair : get_neighbor_list(system)) {
                md::index const i = pair.first;
                md::index const j = pair.second;

                auto const pot = derived().neighbor_pair_potential(system, i, j);
                auto const r = positions[i] - positions[j];

                sum += pot.evaluate_energy(r);
            }

            return sum;
        }

        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view(md::position_attribute);

            for (auto const pair : get_neighbor_list(system)) {
                md::index const i = pair.first;
                md::index const j = pair.second;

                auto const pot = derived().neighbor_pair_potential(system, i, j);
                auto const r = positions[i] - positions[j];

                auto const force = pot.evaluate_force(r);
                forces[i] += force;
                forces[j] -= force;
            }
        }

    private:
        md::neighbor_list& get_neighbor_list(md::system const& system)
        {
            neighbor_list_.update(
                system.view(md::position_attribute),
                derived().neighbor_distance(system)
            );
            return neighbor_list_;
        }

        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

    private:
        md::neighbor_list neighbor_list_;
    };
}

#endif
