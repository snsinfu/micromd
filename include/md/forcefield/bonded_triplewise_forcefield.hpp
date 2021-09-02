// Copyright snsinfu 2019-2021.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_BONDED_TRIPLEWISE_FORCEFIELD_HPP
#define MD_FORCEFIELD_BONDED_TRIPLEWISE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// interactions among selected particle triples.

#include <utility>
#include <vector>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/triple_potfun.hpp"


namespace md
{
    // bonded_triplewise_forcefield implements md::forcefield. It computes
    // interactions between selected triples of particles.
    //
    // This is a CRTP base class. Derived class must define a callback:
    //
    //     auto bonded_triplewise_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) triple.
    //
    template<typename Derived>
    class bonded_triplewise_forcefield : public virtual md::forcefield
    {
    public:
        // add_bonded_triple selects given triple as interacting.
        Derived& add_bonded_triple(md::index i, md::index j, md::index k)
        {
            triples_.push_back({ i, j, k });
            return derived();
        }

        // add_bonded_range selects all adjacent triples in the range [start,end)
        // as interacting.
        Derived& add_bonded_range(md::index start, md::index end)
        {
            for (md::index i = start; i + 2 < end; i++) {
                triples_.push_back({ i, i + 1, i + 2 });
            }
            return derived();
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (index_triple const& triple : triples_) {
                md::index const i = triple.i;
                md::index const j = triple.j;
                md::index const k = triple.k;
                md::vector const rij = positions[i] - positions[j];
                md::vector const rjk = positions[j] - positions[k];

                md::scalar const energy = derived()
                    .bonded_triplewise_potential(system, i, j, k)
                    .evaluate_energy(rij, rjk);

                sum += energy;
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (index_triple const& triple : triples_) {
                md::index const i = triple.i;
                md::index const j = triple.j;
                md::index const k = triple.k;
                md::vector const rij = positions[i] - positions[j];
                md::vector const rjk = positions[j] - positions[k];

                auto const force = derived()
                    .bonded_triplewise_potential(system, i, j, k)
                    .evaluate_force(rij, rjk);

                forces[i] += std::get<0>(force);
                forces[j] += std::get<1>(force);
                forces[k] += std::get<2>(force);
            }
        }

    private:
        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        struct index_triple
        {
            md::index i;
            md::index j;
            md::index k;
        };

        std::vector<index_triple> triples_;
    };

    template<typename PotFun>
    class basic_bonded_triplewise_forcefield
        : public md::bonded_triplewise_forcefield<basic_bonded_triplewise_forcefield<PotFun>>
    {
    public:
        explicit basic_bonded_triplewise_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto bonded_triplewise_potential(
            md::system const& system, md::index i, md::index j, md::index k
        ) const
        {
            return potfun_(system, i, j, k);
        }

    private:
        PotFun potfun_;
    };

    // make_bonded_triplewise_forcefield implements
    // md::bonded_triplewise_forcefield with given potential object or lambda
    // returning a potential object.
    template<typename P>
    auto make_bonded_triplewise_forcefield(P pot)
    {
        auto potfun = detail::make_triple_potential_factory(pot);
        using potfun_type = decltype(potfun);
        return md::basic_bonded_triplewise_forcefield<potfun_type>{potfun};
    }
}

#endif
