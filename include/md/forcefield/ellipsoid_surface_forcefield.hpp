// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_ELLIPSOID_SURFACE_FORCEFIELD_HPP
#define MD_FORCEFIELD_ELLIPSOID_SURFACE_FORCEFIELD_HPP

// This module provides a forcefield template to force particles near surface of
// an ellipsoid.
//
// struct ellipsoid
// class  ellipsoid_surface_forcefield

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "../potential/constant_potential.hpp"


namespace md
{
    // ellipsoid is a triaxial ellipsoid.
    struct ellipsoid
    {
        // Center of the ellipsoid. Defaults to the origin.
        md::point center;

        // Semiaxis length along the x axis. Defaults to 1.
        md::scalar semiaxis_x = 1;

        // Semiaxis length along the y axis. Defaults to 1.
        md::scalar semiaxis_y = 1;

        // Semiaxis length along the z axis. Defaults to 1.
        md::scalar semiaxis_z = 1;
    };

    namespace detail
    {
        struct ellipsoid_eval
        {
            bool undefined;
            md::vector delta;
            md::vector strain;
            md::scalar implicit;
        };

        inline detail::ellipsoid_eval evaluate_point(md::ellipsoid e, md::point pt)
        {
            md::vector const quadform = {
                1 / (e.semiaxis_x * e.semiaxis_x),
                1 / (e.semiaxis_y * e.semiaxis_y),
                1 / (e.semiaxis_z * e.semiaxis_z)
            };

            md::vector const radial = pt - e.center;
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
    template<typename Derived>
    class ellipsoid_surface_forcefield : public virtual md::forcefield
    {
    public:
        // set_ellipsoid changes the ellipsoid to given one.
        void set_ellipsoid(md::ellipsoid ellipsoid)
        {
            ellipsoid_ = ellipsoid;
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (md::index i = 0; i < system.particle_count(); i++) {
                detail::ellipsoid_eval const ev = detail::evaluate_point(ellipsoid_, positions[i]);

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
            md::array_view<md::point const> positions = system.view_positions();

            for (md::index i = 0; i < system.particle_count(); i++) {
                detail::ellipsoid_eval const ev = detail::evaluate_point(ellipsoid_, positions[i]);

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
            }
        }

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

    private:

        // derived returns a reference to this as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        md::ellipsoid ellipsoid_;
    };
}

#endif
