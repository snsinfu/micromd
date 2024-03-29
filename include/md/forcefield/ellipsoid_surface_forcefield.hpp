// Copyright snsinfu 2018-2021.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_ELLIPSOID_SURFACE_FORCEFIELD_HPP
#define MD_FORCEFIELD_ELLIPSOID_SURFACE_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// field force acting on particles near an ellipsoidal surface.

#include <functional>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "../potential/constant_potential.hpp"

#include "detail/field_potfun.hpp"


namespace md
{
    // ellipsoid is a triaxial ellipsoid.
    struct ellipsoid
    {
        // Center of the ellipsoid. Defaults to the origin.
        md::point center;

        // Semiaxis lengths along the coordinate axes. Defaults to 1.
        md::scalar semiaxis_x = 1;
        md::scalar semiaxis_y = 1;
        md::scalar semiaxis_z = 1;

        // implicit returns the implicit ellipsoid function evaluated at pt:
        // f(x,y,z) = (x/a)^2 + (y/b)^2 + (z/c)^2 - 1.
        md::scalar implicit(md::point pt) const
        {
            md::vector const quadform = {
                1 / (semiaxis_x * semiaxis_x),
                1 / (semiaxis_y * semiaxis_y),
                1 / (semiaxis_z * semiaxis_z)
            };
            md::vector const r = pt - center;

            return quadform.hadamard(r).dot(r) - 1;
        }
    };

    namespace detail
    {
        // ellipsoid_eval holds mathematical quantities required to compute the
        // gradient of surface potential on a given point.
        struct ellipsoid_eval
        {
            // True if the point is at the center. Other members are undefined
            // if this flag is true.
            bool undefined;

            // Vector pointing to the point from an estimated nearest point on
            // the ellipsoid.
            md::vector delta;

            // Some vector needed to transform a gradient.
            md::vector strain;

            // The value of the implicit function defining the ellipsoid.
            md::scalar implicit;
        };

        inline detail::ellipsoid_eval evaluate_point(md::ellipsoid ellip, md::point pt)
        {
            md::vector const quadform = {
                1 / (ellip.semiaxis_x * ellip.semiaxis_x),
                1 / (ellip.semiaxis_y * ellip.semiaxis_y),
                1 / (ellip.semiaxis_z * ellip.semiaxis_z)
            };

            md::vector const radial = pt - ellip.center;
            md::vector const dual = quadform.hadamard(radial);

            if (dual.squared_norm() == 0) {
                detail::ellipsoid_eval ev;
                ev.undefined = true;
                return ev;
            }

            md::scalar const implicit = dual.dot(radial) - 1;
            md::scalar const scale = implicit / (2 * dual.squared_norm());

            detail::ellipsoid_eval ev;

            ev.undefined = false;
            ev.delta = scale * dual;
            ev.strain = scale * quadform;
            ev.implicit = implicit;

            return ev;
        }
    }

    // ellipsoid_surface_forcefield computes short-range field interaction of
    // particles and an ellipsoidal surface. It uses an approximation and so
    // computed interactions are inaccurate near the center.
    //
    // This is a CRTP base class. Callbacks are:
    //
    //     auto ellipsoid_inward_potential(
    //         md::system const& system,
    //         md::index i
    //     )
    //     Returns the potential object for a particle inside ellipsoid. It
    //     defaults to a zero potential if not defined.
    //
    //     auto ellipsoid_outward_potential(
    //         md::system const& system,
    //         md::index i
    //     )
    //     Returns the potential object for a particle outside ellipsoid. It
    //     defaults to a zero potential if not defined.
    //
    template<typename Derived>
    class ellipsoid_surface_forcefield : public virtual md::forcefield
    {
    public:
        struct statistics
        {
            // Sum of the normal reaction force acting on the surface calculated
            // in the previous call of compute_force().
            md::scalar reaction_force = 0;

            // Sum of axis-projected reaction force acting on the surface.
            md::vector axial_reaction;
        };
        statistics stats;

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::ellipsoid const ellipsoid = derived().ellipsoid(system);
            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (md::index i = 0; i < system.particle_count(); i++) {
                detail::ellipsoid_eval const ev = detail::evaluate_point(ellipsoid, positions[i]);

                if (ev.undefined) {
                    continue;
                }

                if (ev.implicit < 0) {
                    sum += derived()
                        .ellipsoid_inward_potential(system, i)
                        .evaluate_energy(ev.delta);
                } else {
                    sum += derived()
                        .ellipsoid_outward_potential(system, i)
                        .evaluate_energy(ev.delta);
                }
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::ellipsoid const ellipsoid = derived().ellipsoid(system);
            md::array_view<md::point const> positions = system.view_positions();

            stats.reaction_force = 0;
            stats.axial_reaction = {};

            for (md::index i = 0; i < system.particle_count(); i++) {
                detail::ellipsoid_eval const ev = detail::evaluate_point(ellipsoid, positions[i]);

                if (ev.undefined) {
                    continue;
                }

                md::vector basic_force;

                if (ev.implicit < 0) {
                    basic_force = derived()
                        .ellipsoid_inward_potential(system, i)
                        .evaluate_force(ev.delta);
                } else {
                    basic_force = derived()
                        .ellipsoid_outward_potential(system, i)
                        .evaluate_force(ev.delta);
                }

                md::vector const iso_force = basic_force.project(ev.delta);
                md::vector const aniso_force = basic_force - iso_force;
                md::vector const strain_force = (aniso_force - iso_force).hadamard(ev.strain);
                md::vector const force = iso_force + strain_force;

                forces[i] += force;

                // Collect statistics.
                md::vector const normal =
                    (ev.implicit < 0 ? -1 : 1) * ev.delta.normalize();

                stats.reaction_force -= force.dot(normal);

                stats.axial_reaction.x -= (normal.x > 0 ? force.x : -force.x);
                stats.axial_reaction.y -= (normal.y > 0 ? force.y : -force.y);
                stats.axial_reaction.z -= (normal.z > 0 ? force.z : -force.z);
            }
        }

        //
        // CRTP default implementations
        //

        // ellipsoid_inward_potential by default returns a zero potential.
        md::constant_potential ellipsoid_inward_potential(md::system const&, md::index) const
        {
            return md::constant_potential{0};
        }

        // ellipsoid_outward_potential by default returns a zero potential.
        md::constant_potential ellipsoid_outward_potential(md::system const&, md::index) const
        {
            return md::constant_potential{0};
        }

        // ellipsoid by default returns a unit sphere centered at the origin.
        md::ellipsoid ellipsoid(md::system const&) const
        {
            return {};
        }

    private:
        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }
    };


    template<typename Derived>
    class basic_ellipsoid_surface_forcefield
        : public md::ellipsoid_surface_forcefield<Derived>
    {
    public:
        md::ellipsoid ellipsoid(md::system const&) const
        {
            return ellipsoid_callback_();
        }

        Derived& set_ellipsoid(md::ellipsoid el)
        {
            return set_ellipsoid([=] { return el; });
        }

        Derived& set_ellipsoid(std::function<md::ellipsoid()> el_cb)
        {
            ellipsoid_callback_ = el_cb;
            return static_cast<Derived&>(*this);
        }

    private:
        std::function<md::ellipsoid()> ellipsoid_callback_ = [] { return md::ellipsoid {}; };
    };


    template<typename PotFun>
    class basic_ellipsoid_inward_forcefield_impl
        : public md::basic_ellipsoid_surface_forcefield<
            basic_ellipsoid_inward_forcefield_impl<PotFun>
        >
    {
    public:
        explicit basic_ellipsoid_inward_forcefield_impl(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto ellipsoid_inward_potential(md::system const& system, md::index i) const
        {
            return potfun_(system, i);
        }

    private:
        PotFun potfun_;
    };


    template<typename PotFun>
    class basic_ellipsoid_outward_forcefield_impl
        : public md::basic_ellipsoid_surface_forcefield<basic_ellipsoid_outward_forcefield_impl<PotFun>>
    {
    public:
        explicit basic_ellipsoid_outward_forcefield_impl(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto ellipsoid_outward_potential(md::system const& system, md::index i) const
        {
            return potfun_(system, i);
        }

    private:
        PotFun potfun_;
    };


    // make_ellipsoid_inward_forcefield implements
    // md::ellipsoid_surface_forcefield with given potential object or lambda
    // returning a potential object.
    template<typename P>
    auto make_ellipsoid_inward_forcefield(P pot)
    {
        auto potfun = detail::make_field_potential_factory(pot);
        using potfun_type = decltype(potfun);
        return md::basic_ellipsoid_inward_forcefield_impl<potfun_type>{potfun};
    }


    // make_ellipsoid_outward_forcefield implements
    // md::ellipsoid_surface_forcefield with given potential object or lambda
    // returning a potential object.
    template<typename P>
    auto make_ellipsoid_outward_forcefield(P pot)
    {
        auto potfun = detail::make_field_potential_factory(pot);
        using potfun_type = decltype(potfun);
        return md::basic_ellipsoid_outward_forcefield_impl<potfun_type>{potfun};
    }
}

#endif
