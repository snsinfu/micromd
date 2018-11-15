// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_NEIGHBOR_PAIR_FORCEFIELD_HPP
#define MD_FORCEFIELD_NEIGHBOR_PAIR_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// short-range interactions between particles in a cutoff distance.

#include <type_traits>
#include <utility>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/neighbor_list.hpp"


namespace md
{
    // neighbor_pair_forcefield implements md::forcefield. It computes short-
    // range interactions between every pair of particles that are close within
    // a given cutoff distance.
    //
    // This is a CRTP base class. Derived class must define two callbacks:
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
    template<typename Derived>
    class neighbor_pair_forcefield : public virtual md::forcefield
    {
    public:
        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();
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

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

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
        // get_neighbor_list returns a reference to the up-to-date neighbor list
        // for the system.
        md::neighbor_list& get_neighbor_list(md::system const& system)
        {
            neighbor_list_.update(
                system.view(md::position_attribute),
                derived().neighbor_distance(system)
            );
            return neighbor_list_;
        }

        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        md::neighbor_list neighbor_list_;
    };

    // Undocumented experimental features:

    namespace detail
    {
        template<typename PF>
        class basic_neighbor_pair_forcefield
            : public md::neighbor_pair_forcefield<basic_neighbor_pair_forcefield<PF>>
        {
        public:
            basic_neighbor_pair_forcefield(md::scalar dcut, PF potfun)
                : dcut_{dcut}, potfun_{potfun}
            {
            }

            inline md::scalar neighbor_distance(md::system const&) const
            {
                return dcut_;
            }

            inline auto neighbor_pair_potential(md::system const&, md::index i, md::index j) const
            {
                return potfun_(i, j);
            }

        private:
            md::scalar dcut_;
            PF potfun_;
        };
    }

    namespace detail
    {
        template<typename P>
        struct pair_potential_factory
        {
            P potential;

            P operator()(md::index, md::index) const
            {
                return potential;
            }
        };

        template<
            typename PF,
            typename = decltype(std::declval<PF>()(md::index{}, md::index{}))
        >
        PF make_pair_potfun(PF potfun)
        {
            return potfun;
        }

        template<typename P, typename... Dummy>
        pair_potential_factory<P> make_pair_potfun(P pot, Dummy...)
        {
            return pair_potential_factory<P>{pot};
        }
    }

    template<typename P>
    void force_neighbor_pairs(md::system& system, md::scalar dcut, P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        using potfun_t = decltype(potfun);

        system.add_forcefield(
            std::make_shared<detail::basic_neighbor_pair_forcefield<potfun_t>>(dcut, potfun)
        );
    }
}

#endif
