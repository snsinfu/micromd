// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_POINT_SOURCE_FORCEFIELD_HPP
#define MD_FORCEFIELD_POINT_SOURCE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes field
// force from a point source.

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"


namespace md
{
    // point_source_forcefield computes field interaction of particles and a
    // fixed point source.
    //
    // This is a CRTP base class. Derived class must implement a callback:
    //
    //     auto point_source_potential(
    //         md::system const& system,
    //         md::index i
    //     )
    //     Returns the potential object for a particle.
    //
    template<typename Derived>
    class point_source_forcefield : public virtual md::forcefield
    {
    public:
        // set_point_source sets a source point.
        Derived& set_point_source(md::point pt)
        {
            source_ = pt;
            return derived();
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = positions[i] - source_;

                sum += derived()
                    .point_source_potential(system, i)
                    .evaluate_energy(r);
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = positions[i] - source_;

                forces[i] += derived()
                    .point_source_potential(system, i)
                    .evaluate_force(r);
            }
        }

    private:
        // derived returns a reference to this as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        md::point source_;
    };
}

#endif
