// Copyright snsinfu 2018-2021.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_BRUTEFORCE_PAIRWISE_FORCEFIELD_HPP
#define MD_FORCEFIELD_BRUTEFORCE_PAIRWISE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// interactions between all particle pairs.

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/pair_potfun.hpp"


namespace md
{
    // bruteforce_pairwise_forcefield implements md::forcefield. It computes
    // interactions between every pair of particles.
    //
    // This is a CRTP base class. Derived class must define a callback:
    //
    //     auto bruteforce_pairwise_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) pair.
    //
    template<typename Derived>
    class bruteforce_pairwise_forcefield : public virtual md::forcefield
    {
    public:
        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();
            md::scalar sum = 0;

            for (md::index j = 0; j < positions.size(); j++) {
                for (md::index i = 0; i < j; i++) {
                    auto const pot = derived().bruteforce_pairwise_potential(system, i, j);
                    auto const r = positions[i] - positions[j];

                    sum += pot.evaluate_energy(r);
                }
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (md::index j = 0; j < positions.size(); j++) {
                for (md::index i = 0; i < j; i++) {
                    auto const pot = derived().bruteforce_pairwise_potential(system, i, j);
                    auto const r = positions[i] - positions[j];

                    auto const force = pot.evaluate_force(r);
                    forces[i] += force;
                    forces[j] -= force;
                }
            }
        }

    private:
        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }
    };

    template<typename PotFun>
    class basic_bruteforce_pairwise_forcefield
        : public md::bruteforce_pairwise_forcefield<basic_bruteforce_pairwise_forcefield<PotFun>>
    {
    public:
        explicit basic_bruteforce_pairwise_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto bruteforce_pairwise_potential(md::system const& system, md::index i, md::index j) const
        {
            return potfun_(system, i, j);
        }

    private:
        PotFun potfun_;
    };

    // make_bruteforce_pairwise_forcefield implements
    // md::bruteforce_pairwise_forcefield with given potential object or lambda
    // returning a potential object.
    template<typename P>
    auto make_bruteforce_pairwise_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potential_factory(pot);
        using potfun_type = decltype(potfun);
        return md::basic_bruteforce_pairwise_forcefield<potfun_type>{potfun};
    }
}

#endif
