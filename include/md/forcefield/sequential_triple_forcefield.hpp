// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_SEQUENTIAL_TRIPLE_FORCEFIELD_HPP
#define MD_FORCEFIELD_SEQUENTIAL_TRIPLE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// interactions among sequentially numbered particle triplets.

#include <tuple>
#include <utility>
#include <vector>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/triple_potfun.hpp"


namespace md
{
    // sequential_triple_forcefield implements md::forcefield. It computes
    // interactions among every consecutive triple of particles in given
    // segments.
    //
    // This is a CRTP base class. Derived class must define a callback:
    //
    //     auto sequential_triple_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j,
    //         md::index k
    //     )
    //     Returns the potential object for (i,j,k) pair.
    //
    template<typename Derived>
    class sequential_triple_forcefield : public virtual md::forcefield
    {
    public:
        // add_segment marks all consecutive particle triplets in a segment as
        // interacting. The range [first,last] is inclusive.
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

            for (std::pair<md::index, md::index> const segment : segments_) {
                for (md::index i = segment.first; i < segment.second - 1; i++) {
                    md::index const j = i + 1;
                    md::index const k = i + 2;
                    md::vector const rij = positions[i] - positions[j];
                    md::vector const rjk = positions[j] - positions[k];

                    md::scalar const energy = derived()
                        .sequential_triple_potential(system, i, j, k)
                        .evaluate_energy(rij, rjk);

                    sum += energy;
                }
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (std::pair<md::index, md::index> const segment : segments_) {
                for (md::index i = segment.first; i < segment.second - 1; i++) {
                    md::index const j = i + 1;
                    md::index const k = i + 2;
                    md::vector const rij = positions[i] - positions[j];
                    md::vector const rjk = positions[j] - positions[k];

                    auto const force = derived()
                        .sequential_triple_potential(system, i, j, k)
                        .evaluate_force(rij, rjk);

                    forces[i] += std::get<0>(force);
                    forces[j] += std::get<1>(force);
                    forces[k] += std::get<2>(force);
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
    class basic_sequential_triple_forcefield
        : public md::sequential_triple_forcefield<basic_sequential_triple_forcefield<PotFun>>
    {
    public:
        explicit basic_sequential_triple_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto sequential_triple_potential(
            md::system const&, md::index i, md::index j, md::index k
        ) const
        {
            return potfun_(i, j, k);
        }

    private:
        PotFun potfun_;
    };

    // make_sequential_triple_forcefield implements
    // md::sequential_triple_forcefield with given potential object or lambda
    // returning a potential object.
    template<typename P>
    auto make_sequential_triple_forcefield(P pot)
    {
        auto potfun = detail::make_triple_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_sequential_triple_forcefield<potfun_type>{potfun};
    }
}

#endif
