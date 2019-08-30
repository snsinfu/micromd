// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_NEIGHBOR_PAIR_FORCEFIELD_V2_HPP
#define MD_FORCEFIELD_NEIGHBOR_PAIR_FORCEFIELD_V2_HPP

// This module provides a template forcefield implementation that quickly
// computes short-range pairwise interactions in open and periodic systems.

#include <type_traits>
#include <utility>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/neighbor_list_v2.hpp"
#include "detail/pair_potfun.hpp"


namespace md
{
    // neighbor_pair_forcefield_v2 implements md::forcefield. It computes short-
    // range interactions between every pair of particles that are close within
    // a given cutoff distance.
    //
    // This is a CRTP base class. Derived class must define three callbacks:
    //
    //     Box box(md::system const& system)
    //     Returns the box description.
    //
    //     md::scalar neighbor_distance(md::system const& system)
    //     Returns the cutoff distance.
    //
    //     auto neighbor_pair_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) pair.
    //
    // TODO: Reconsider the name of the box() callback.
    //
    template<typename Box, typename Derived>
    class neighbor_pair_forcefield_v2 : public virtual md::forcefield
    {
    public:
        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            Box const box = derived().box(system);
            md::array_view<md::point const> positions = system.view_positions();
            md::scalar sum = 0;

            for (auto const pair : get_neighbor_list(system)) {
                md::index const i = pair.first;
                md::index const j = pair.second;

                auto const pot = derived().neighbor_pair_potential(system, i, j);
                auto const r = box.shortest_displacement(positions[i], positions[j]);

                sum += pot.evaluate_energy(r);
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            Box const box = derived().box(system);
            md::array_view<md::point const> positions = system.view_positions();

            for (auto const pair : get_neighbor_list(system)) {
                md::index const i = pair.first;
                md::index const j = pair.second;

                auto const pot = derived().neighbor_pair_potential(system, i, j);
                auto const r = box.shortest_displacement(positions[i], positions[j]);

                auto const force = pot.evaluate_force(r);
                forces[i] += force;
                forces[j] -= force;
            }
        }

    private:
        // get_neighbor_list returns a reference to the up-to-date neighbor list
        // for the system.
        md::neighbor_list_v2<Box> const& get_neighbor_list(md::system const& system)
        {
            neighbor_list_.update(
                system.view_positions(),
                derived().neighbor_distance(system),
                derived().box(system)
            );
            return neighbor_list_;
        }

        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        md::neighbor_list_v2<Box> neighbor_list_;
    };

    // Implements md::neighbor_pair_forcefield_v2. See the factory function
    // make_neighbor_pair_forcefield_v2.
    template<typename Box, typename PotFun>
    class basic_neighbor_pair_forcefield_v2
        : public md::neighbor_pair_forcefield_v2<
            Box, basic_neighbor_pair_forcefield_v2<Box, PotFun>
        >
    {
    public:
        explicit basic_neighbor_pair_forcefield_v2(PotFun potfun)
            : potfun_{potfun}
        {
        }

        basic_neighbor_pair_forcefield_v2& set_box(Box box)
        {
            box_ = box;
            return *this;
        }

        basic_neighbor_pair_forcefield_v2& set_neighbor_distance(md::scalar ndist)
        {
            ndist_ = ndist;
            return *this;
        }

        Box box(md::system const&) const
        {
            return box_;
        }

        md::scalar neighbor_distance(md::system const&) const
        {
            return ndist_;
        }

        auto neighbor_pair_potential(md::system const&, md::index i, md::index j) const
        {
            return potfun_(i, j);
        }

    private:
        Box box_;
        PotFun potfun_;
        md::scalar ndist_ = 0;
    };

    // make_neighbor_pair_forcefield_v2 implements md::neighbor_pair_forcefield_v2
    // with given potential object or lambda returning a potential object.
    template<typename Box = md::open_box, typename P>
    auto make_neighbor_pair_forcefield_v2(P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_neighbor_pair_forcefield_v2<Box, potfun_type>{potfun};
    }
}

#endif
