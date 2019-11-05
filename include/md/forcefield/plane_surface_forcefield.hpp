// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_PLANE_SURFACE_FORCEFIELD_HPP
#define MD_FORCEFIELD_PLANE_SURFACE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// field force acting on particles near a plane.

#include <cmath>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "../potential/constant_potential.hpp"

#include "detail/field_potfun.hpp"


namespace md
{
    // plane_surface_forcefield computes field interaction of particles and a
    // plane.
    //
    // This is a CRTP base class. Callbacks are:
    //
    //     auto plane_inward_potential(
    //         md::system const& system,
    //         md::index i
    //     )
    //     Returns the potential object for a particle on the negative side of
    //     the plane. It defaults to a zero potential if not defined.
    //
    //     auto plane_outward_potential(
    //         md::system const& system,
    //         md::index i
    //     )
    //     Returns the potential object for a particle on the positive side of
    //     the plane. It defaults to a zero potential if not defined.
    //
    template<typename Derived>
    class plane_surface_forcefield : public virtual md::forcefield
    {
    public:
        struct statistics
        {
            // Sum of the normal reaction force acting on the surface calculated
            // in the previous call of compute_force().
            md::scalar reaction_force = 0;
        };
        statistics stats;

        // set_plane changes the plane. Default is the xy plane.
        Derived& set_plane(md::point reference, md::vector normal)
        {
            reference_ = reference;
            normal_ = normal;
            return derived();
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();
            md::scalar sum = 0;

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = (positions[i] - reference_).project(normal_);

                if (r.dot(normal_) < 0) {
                    auto const pot = derived().plane_inward_potential(system, i);
                    sum += pot.evaluate_energy(r);
                } else {
                    auto const pot = derived().plane_outward_potential(system, i);
                    sum += pot.evaluate_energy(r);
                }
            }
            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = (positions[i] - reference_).project(normal_);

                if (r.dot(normal_) < 0) {
                    auto const pot = derived().plane_inward_potential(system, i);
                    forces[i] += pot.evaluate_force(r);
                } else {
                    auto const pot = derived().plane_outward_potential(system, i);
                    forces[i] += pot.evaluate_force(r);
                }
            }
        }

        // plane_inward_potential by default returns a zero potential.
        md::constant_potential plane_inward_potential(md::system const&, md::index) const
        {
            return md::constant_potential{0};
        }

        // plane_outward_potential by default returns a zero potential.
        md::constant_potential plane_outward_potential(md::system const&, md::index) const
        {
            return md::constant_potential{0};
        }

    private:
        // derived returns a reference to this as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        md::point reference_ = {0, 0, 0};
        md::vector normal_ = {0, 0, 1};
    };

    template<typename PotFun>
    class basic_plane_inward_forcefield
        : public md::plane_surface_forcefield<basic_plane_inward_forcefield<PotFun>>
    {
    public:
        explicit basic_plane_inward_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto plane_inward_potential(md::system const&, md::index i) const
        {
            return potfun_(i);
        }

    private:
        PotFun potfun_;
    };

    template<typename PotFun>
    class basic_plane_outward_forcefield
        : public md::plane_surface_forcefield<basic_plane_outward_forcefield<PotFun>>
    {
    public:
        explicit basic_plane_outward_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto plane_outward_potential(md::system const&, md::index i) const
        {
            return potfun_(i);
        }

    private:
        PotFun potfun_;
    };

    // make_plane_inward_forcefield implements md::plane_surface_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename P>
    auto make_plane_inward_forcefield(P pot)
    {
        auto potfun = detail::make_field_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_plane_inward_forcefield<potfun_type>{potfun};
    }

    // make_plane_outward_forcefield implements md::plane_surface_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename P>
    auto make_plane_outward_forcefield(P pot)
    {
        auto potfun = detail::make_field_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_plane_outward_forcefield<potfun_type>{potfun};
    }
}

#endif
