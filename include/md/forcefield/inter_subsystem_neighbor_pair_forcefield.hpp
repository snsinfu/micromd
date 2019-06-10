// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_INTER_SUBSYSTEM_NEIGHBOR_PAIR_FORCEFIELD_HPP
#define MD_FORCEFIELD_INTER_SUBSYSTEM_NEIGHBOR_PAIR_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// short-range interactions between neighboring particles in two subsystems.

#include <type_traits>
#include <utility>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/inter_subsystem_neighbor_list.hpp"
#include "detail/pair_potfun.hpp"


namespace md
{
    // inter_subsystem_neighbor_pair_forcefield implements md::forcefield. It
    // computes short-range interactions between every pair of particles in two
    // subsystems that are closer than a given cutoff distance.
    //
    // This is a CRTP base class. Derived class must define two callbacks:
    //
    //     md::scalar inter_subsystem_neighbor_distance(md::system const& system)
    //     Returns the cutoff distance.
    //
    //     auto inter_subsystem_neighbor_pair_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) pair.
    //
    template<typename Derived>
    class inter_subsystem_neighbor_pair_forcefield : public virtual md::forcefield
    {
    public:
        // set_subsystems sets two subsystems to operate on.
        Derived& set_subsystems(
            md::array_view<md::index const> subsystem1,
            md::array_view<md::index const> subsystem2
        )
        {
            neighbor_list_.set_subsystems(subsystem1, subsystem2);
            return derived();
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();
            md::scalar sum = 0;

            for (auto const pair : get_neighbor_list(system)) {
                md::index const i = pair.first;
                md::index const j = pair.second;

                auto const pot = derived().inter_subsystem_neighbor_pair_potential(
                    system, i, j
                );
                auto const r = positions[i] - positions[j];

                sum += pot.evaluate_energy(r);
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (auto const pair : get_neighbor_list(system)) {
                md::index const i = pair.first;
                md::index const j = pair.second;

                auto const pot = derived().inter_subsystem_neighbor_pair_potential(
                    system, i, j
                );
                auto const r = positions[i] - positions[j];

                auto const force = pot.evaluate_force(r);
                forces[i] += force;
                forces[j] -= force;
            }
        }

        // get_neighbor_list returns a reference to the up-to-date neighbor list
        // for the subsystem.
        md::inter_subsystem_neighbor_list const& get_neighbor_list(md::system const& system)
        {
            neighbor_list_.update(
                system.view_positions(),
                derived().inter_subsystem_neighbor_distance(system)
            );
            return neighbor_list_;
        }

    private:
        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        md::inter_subsystem_neighbor_list neighbor_list_;
        std::vector<md::index> key_indices_;
    };

    template<typename PotFun>
    class basic_inter_subsystem_neighbor_pair_forcefield
        : public md::inter_subsystem_neighbor_pair_forcefield<
            basic_inter_subsystem_neighbor_pair_forcefield<PotFun>
        >
    {
    public:
        explicit basic_inter_subsystem_neighbor_pair_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        basic_inter_subsystem_neighbor_pair_forcefield& set_neighbor_distance(md::scalar ndist)
        {
            ndist_ = ndist;
            return *this;
        }

        md::scalar inter_subsystem_neighbor_distance(md::system const&) const
        {
            return ndist_;
        }

        auto inter_subsystem_neighbor_pair_potential(
            md::system const&, md::index i, md::index j
        ) const
        {
            return potfun_(i, j);
        }

    private:
        PotFun potfun_;
        md::scalar ndist_;
    };

    // make_inter_subsystem_neighbor_pair_forcefield implements
    // md::inter_subsystem_neighbor_pair_forcefield with given potential object
    // or lambda returning a potential object.
    template<typename P>
    auto make_inter_subsystem_neighbor_pair_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_inter_subsystem_neighbor_pair_forcefield<potfun_type>{potfun};
    }
}

#endif
