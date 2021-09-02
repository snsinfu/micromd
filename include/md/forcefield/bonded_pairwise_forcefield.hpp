// Copyright snsinfu 2019-2021.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_BONDED_PAIRWISE_FORCEFIELD_HPP
#define MD_FORCEFIELD_BONDED_PAIRWISE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// interactions between selected particle pairs.

#include <utility>
#include <vector>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/pair_potfun.hpp"


namespace md
{
    // bonded_pairwise_forcefield implements md::forcefield. It computes
    // interactions between selected pairs of particles.
    //
    // This is a CRTP base class. Derived class must define a callback:
    //
    //     auto bonded_pairwise_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) pair.
    //
    template<typename Derived>
    class bonded_pairwise_forcefield : public virtual md::forcefield
    {
    public:
        // add_bonded_pair selects given pair as interacting.
        Derived& add_bonded_pair(md::index i, md::index j)
        {
            pairs_.emplace_back(i, j);
            return derived();
        }

        // add_bonded_range selects all adjacent pairs in the range [start,end)
        // as interacting.
        Derived& add_bonded_range(md::index start, md::index end)
        {
            for (md::index i = start; i + 1 < end; i++) {
                pairs_.emplace_back(i, i + 1);
            }
            return derived();
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (std::pair<md::index, md::index> pair : pairs_) {
                md::index const i = pair.first;
                md::index const j = pair.second;
                md::vector const r = positions[i] - positions[j];

                auto const pot = derived().bonded_pairwise_potential(system, i, j);
                md::scalar const energy = pot.evaluate_energy(r);

                sum += energy;
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (std::pair<md::index, md::index> pair : pairs_) {
                md::index const i = pair.first;
                md::index const j = pair.second;
                md::vector const r = positions[i] - positions[j];

                auto const pot = derived().bonded_pairwise_potential(system, i, j);
                md::vector const force = pot.evaluate_force(r);

                forces[i] += force;
                forces[j] -= force;
            }
        }

    private:
        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        std::vector<std::pair<md::index, md::index>> pairs_;
    };

    template<typename PotFun>
    class basic_bonded_pairwise_forcefield
        : public md::bonded_pairwise_forcefield<basic_bonded_pairwise_forcefield<PotFun>>
    {
    public:
        explicit basic_bonded_pairwise_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto bonded_pairwise_potential(md::system const& system, md::index i, md::index j) const
        {
            return potfun_(system, i, j);
        }

    private:
        PotFun potfun_;
    };

    // make_bonded_pairwise_forcefield implements md::bonded_pairwise_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename P>
    auto make_bonded_pairwise_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potential_factory(pot);
        using potfun_type = decltype(potfun);
        return md::basic_bonded_pairwise_forcefield<potfun_type>{potfun};
    }
}

#endif
