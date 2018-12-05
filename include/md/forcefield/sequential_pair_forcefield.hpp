// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_SEQUENTIAL_PAIR_FORCEFIELD_HPP
#define MD_FORCEFIELD_SEQUENTIAL_PAIR_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// interactions between sequentially numbered particles.

#include <utility>
#include <vector>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/pair_potfun.hpp"


namespace md
{
    // sequential_pair_forcefield implements md::forcefield. It computes
    // interactions between every consecutive pair of particles in given
    // segments.
    //
    // This is a CRTP base class. Derived class must define a callback:
    //
    //     auto sequential_pair_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) pair.
    //
    template<typename Derived>
    class sequential_pair_forcefield : public virtual md::forcefield
    {
    public:
        // add_segment marks all adjacent particles in a segment as interacting.
        // The range [first,last] is inclusive.
        Derived& add_segment(md::index first, md::index last)
        {
            segments_.emplace_back(first, last);
            return derived();
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (std::pair<md::index, md::index> segment : segments_) {
                for (md::index i = segment.first; i < segment.second; i++) {
                    md::index const j = i + 1;
                    md::vector const r = positions[i] - positions[j];

                    auto const pot = derived().sequential_pair_potential(system, i, j);
                    md::scalar const energy = pot.evaluate_energy(r);

                    sum += energy;
                }
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (std::pair<md::index, md::index> segment : segments_) {
                for (md::index i = segment.first; i < segment.second; i++) {
                    md::index const j = i + 1;
                    md::vector const r = positions[i] - positions[j];

                    auto const pot = derived().sequential_pair_potential(system, i, j);
                    md::vector const force = pot.evaluate_force(r);

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

        std::vector<std::pair<md::index, md::index>> segments_;
    };

    template<typename PotFun>
    class basic_sequential_pair_forcefield
        : public md::sequential_pair_forcefield<basic_sequential_pair_forcefield<PotFun>>
    {
    public:
        explicit basic_sequential_pair_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto sequential_pair_potential(md::system const&, md::index i, md::index j) const
        {
            return potfun_(i, j);
        }

    private:
        PotFun potfun_;
    };

    // make_sequential_pair_forcefield implements md::sequential_pair_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename P>
    auto make_sequential_pair_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_sequential_pair_forcefield<potfun_type>{potfun};
    }
}

#endif
