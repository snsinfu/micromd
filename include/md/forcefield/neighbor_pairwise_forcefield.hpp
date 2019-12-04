// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_NEIGHBOR_PAIRWISE_FORCEFIELD_HPP
#define MD_FORCEFIELD_NEIGHBOR_PAIRWISE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that quickly
// computes short-range pairwise interactions in open and periodic systems.

#include <functional>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"
#include "../misc/index_range.hpp"

#include "detail/neighbor_list.hpp"
#include "detail/pair_potfun.hpp"


namespace md
{
    // neighbor_pairwise_forcefield implements md::forcefield. It computes
    // short- range interactions between every pair of particles that are close
    // within a given cutoff distance.
    //
    // This is a CRTP base class. Derived class must define three callbacks:
    //
    //     Box unit_cell(md::system const& system)
    //     Returns the unit cell of the system.
    //
    //     md::scalar neighbor_distance(md::system const& system)
    //     Returns the cutoff distance.
    //
    //     auto neighbor_pairwise_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) pair.
    //
    template<typename Derived, typename Box = md::open_box>
    class neighbor_pairwise_forcefield : public virtual md::forcefield
    {
    public:
        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            Box const box = derived().unit_cell(system);
            md::array_view<md::point const> positions = system.view_positions();
            md::scalar sum = 0;

            for (auto const pair : get_neighbor_list(system)) {
                md::index const i = pair.first;
                md::index const j = pair.second;

                auto const pot = derived().neighbor_pairwise_potential(system, i, j);
                auto const r = box.shortest_displacement(positions[i], positions[j]);

                sum += pot.evaluate_energy(r);
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            Box const box = derived().unit_cell(system);
            md::array_view<md::point const> positions = system.view_positions();

            for (auto const pair : get_neighbor_list(system)) {
                md::index const i = pair.first;
                md::index const j = pair.second;

                auto const pot = derived().neighbor_pairwise_potential(system, i, j);
                auto const r = box.shortest_displacement(positions[i], positions[j]);

                auto const force = pot.evaluate_force(r);
                forces[i] += force;
                forces[j] -= force;
            }
        }

        Box unit_cell(md::system const&) const
        {
            return {};
        }

        template<typename R>
        Derived& set_neighbor_targets(R const& indices)
        {
            // FIXME: Move to basic_neighbor_pairwise_forcefield. Need to re-
            //        design neighbor_list.
            neighbor_list_.set_targets(indices);
            return derived();
        }

    private:
        // get_neighbor_list returns a reference to the up-to-date neighbor list
        // for the system.
        md::neighbor_list<Box> const& get_neighbor_list(md::system const& system)
        {
            neighbor_list_.update(
                system.view_positions(),
                derived().neighbor_distance(system),
                derived().unit_cell(system)
            );
            return neighbor_list_;
        }

        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        md::neighbor_list<Box> neighbor_list_;
    };


    // Intermediate CRTP layer providing basic implementation of the parameter
    // callbacks of `md::neighbor_pairwise_forcefield`.
    template<typename Derived, typename Box = md::open_box>
    class basic_neighbor_pairwise_forcefield
        : public md::neighbor_pairwise_forcefield<Derived, Box>
    {
    public:
        Derived& set_unit_cell(Box box)
        {
            return set_unit_cell([=] { return box; });
        }

        Derived& set_unit_cell(std::function<Box()> box_cb)
        {
            box_callback_ = box_cb;
            return derived();
        }

        Derived& set_neighbor_distance(md::scalar ndist)
        {
            return set_neighbor_distance([=] { return ndist; });
        }

        Derived& set_neighbor_distance(std::function<md::scalar()> ndist_cb)
        {
            ndist_callback_ = ndist_cb;
            return derived();
        }

        Box unit_cell(md::system const&) const
        {
            return box_callback_();
        }

        md::scalar neighbor_distance(md::system const&) const
        {
            return ndist_callback_();
        }

    private:
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

    private:
        std::function<Box()> box_callback_ = [] { return Box(); };
        std::function<md::scalar()> ndist_callback_ = [] { return 1e-6; };
    };


    template<typename PotFun, typename Box>
    class basic_neighbor_pairwise_forcefield_impl
        : public md::basic_neighbor_pairwise_forcefield<
            basic_neighbor_pairwise_forcefield_impl<PotFun, Box>, Box
        >
    {
    public:
        explicit basic_neighbor_pairwise_forcefield_impl(PotFun const& potfun)
            : potfun_{potfun}
        {
        }

        auto neighbor_pairwise_potential(md::system const&, md::index i, md::index j) const
        {
            return potfun_(i, j);
        }

    private:
        PotFun potfun_;
    };


    // make_neighbor_pairwise_forcefield implements md::neighbor_pairwise_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename Box = md::open_box, typename P>
    auto make_neighbor_pairwise_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_neighbor_pairwise_forcefield_impl<potfun_type, Box>{potfun};
    }
}

#endif
