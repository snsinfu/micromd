// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_SPHERE_SURFACE_FORCEFIELD_HPP
#define MD_FORCEFIELD_SPHERE_SURFACE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// field force acting on particles near a spherical surface.

#include <cmath>
#include <functional>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "../potential/constant_potential.hpp"

#include "detail/field_potfun.hpp"


namespace md
{
    // sphere is a sphere in the three-dimensional space.
    struct sphere
    {
        // Center. Defaults to the origin.
        md::point center;

        // Radius. Defaults to 1.
        md::scalar radius = 1;

        // implicit returns the implicit sphere function evaluated at pt:
        // f(x,y,z) = x^2 + y^2 + z^2 - R^2.
        md::scalar implicit(md::point pt) const
        {
            return (pt - center).squared_norm() - radius * radius;
        }
    };

    // sphere_surface_forcefield computes field interaction of particles and a
    // spherical surface.
    //
    // This is a CRTP base class. Callbacks are:
    //
    //     auto sphere_inward_potential(
    //         md::system const& system,
    //         md::index i
    //     )
    //     Returns the potential object for a particle inside sphere. It
    //     defaults to a zero potential if not defined.
    //
    //     auto sphere_outward_potential(
    //         md::system const& system,
    //         md::index i
    //     )
    //     Returns the potential object for a particle outside sphere. It
    //     defaults to a zero potential if not defined.
    //
    template<typename Derived>
    class sphere_surface_forcefield : public virtual md::forcefield
    {
    public:
        struct statistics
        {
            // Sum of the normal reaction force acting on the surface calculated
            // in the previous call of compute_force().
            md::scalar reaction_force = 0;
        };
        statistics stats;

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::sphere const sphere = derived().sphere(system);
            md::point const center = sphere.center;
            md::scalar const radius = sphere.radius;
            md::scalar const radius2 = radius * radius;

            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = positions[i] - center;
                md::scalar const r2 = r.squared_norm();
                md::scalar const r1 = std::sqrt(r2);

                if (r2 == 0) {
                    continue;
                }

                md::scalar const scale = radius / r1;
                md::vector const s = r - scale * r;

                if (r2 < radius2) {
                    auto const pot = derived().sphere_inward_potential(system, i);
                    sum += pot.evaluate_energy(s);
                } else {
                    auto const pot = derived().sphere_outward_potential(system, i);
                    sum += pot.evaluate_energy(s);
                }
            }
            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::sphere const sphere = derived().sphere(system);
            md::point const center = sphere.center;
            md::scalar const radius = sphere.radius;
            md::scalar const radius2 = radius * radius;

            md::array_view<md::point const> positions = system.view_positions();

            stats.reaction_force = 0;

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = positions[i] - center;
                md::scalar const r2 = r.squared_norm();
                md::scalar const r1 = std::sqrt(r2);

                if (r2 == 0) {
                    continue;
                }

                md::scalar const scale = radius / r1;
                md::vector const s = r - scale * r;

                md::vector force;

                if (r2 < radius2) {
                    auto const pot = derived().sphere_inward_potential(system, i);
                    force = pot.evaluate_force(s);
                } else {
                    auto const pot = derived().sphere_outward_potential(system, i);
                    force = pot.evaluate_force(s);
                }

                md::vector const aniso = scale * (force.project(r) - force);

                forces[i] += force + aniso;
                stats.reaction_force -= force.dot(r) / r.norm();
            }
        }

        //
        // CRTP default implementations
        //

        // sphere_inward_potential by default returns a zero potential.
        md::constant_potential sphere_inward_potential(md::system const&, md::index) const
        {
            return md::constant_potential{0};
        }

        // sphere_outward_potential by default returns a zero potential.
        md::constant_potential sphere_outward_potential(md::system const&, md::index) const
        {
            return md::constant_potential{0};
        }

        // sphere_surface by default returns a unit sphere centered at the origin.
        md::sphere sphere(md::system const&) const
        {
            return {};
        }

    private:
        // derived returns a reference to this as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }
    };


    template<typename Derived>
    class basic_sphere_surface_forcefield : public md::sphere_surface_forcefield<Derived>
    {
    public:
        md::sphere sphere(md::system const&) const
        {
            return sphere_callback_();
        }

        Derived& set_sphere(md::sphere sphere)
        {
            return set_sphere([=] { return sphere; });
        }

        Derived& set_sphere(std::function<md::sphere()> sphere_cb)
        {
            sphere_callback_ = sphere_cb;
            return static_cast<Derived&>(*this);
        }

    private:
        std::function<md::sphere()> sphere_callback_ = [] { return md::sphere{}; };
    };



    template<typename PotFun>
    class basic_sphere_inward_forcefield_impl
        : public md::basic_sphere_surface_forcefield<basic_sphere_inward_forcefield_impl<PotFun>>
    {
    public:
        explicit basic_sphere_inward_forcefield_impl(PotFun const& potfun)
            : potfun_{potfun}
        {
        }

        auto sphere_inward_potential(md::system const&, md::index i) const
        {
            return potfun_(i);
        }

    private:
        PotFun potfun_;
    };


    template<typename PotFun>
    class basic_sphere_outward_forcefield_impl
        : public md::basic_sphere_surface_forcefield<basic_sphere_outward_forcefield_impl<PotFun>>
    {
    public:
        explicit basic_sphere_outward_forcefield_impl(PotFun const& potfun)
            : potfun_{potfun}
        {
        }

        auto sphere_outward_potential(md::system const&, md::index i) const
        {
            return potfun_(i);
        }

    private:
        PotFun potfun_;
    };


    // make_sphere_inward_forcefield implements md::sphere_surface_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename P>
    auto make_sphere_inward_forcefield(P pot)
    {
        auto potfun = detail::make_field_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_sphere_inward_forcefield_impl<potfun_type>{potfun};
    }


    // make_sphere_outward_forcefield implements md::sphere_surface_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename P>
    auto make_sphere_outward_forcefield(P pot)
    {
        auto potfun = detail::make_field_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_sphere_outward_forcefield_impl<potfun_type>{potfun};
    }
}

#endif
