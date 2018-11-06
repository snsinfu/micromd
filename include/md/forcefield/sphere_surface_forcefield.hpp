// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_SPHERE_SURFACE_FORCEFIELD_HPP
#define MD_FORCEFIELD_SPHERE_SURFACE_FORCEFIELD_HPP

#include <type_traits>
#include <utility>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "../potential/constant_potential.hpp"


namespace md
{
    struct sphere
    {
        md::point center;
        md::scalar radius = 1;
    };

    template<typename Derived>
    class sphere_surface_forcefield : public virtual md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const& system) override
        {
            md::sphere const sphere = derived().sphere_surface(system);
            md::scalar const radius2 = sphere.radius * sphere.radius;
            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = positions[i] - sphere.center;
                md::scalar const r2 = r.squared_norm();

                if (r2 < radius2) {
                    auto const pot = derived().sphere_inward_potential(system, i);
                    sum += pot.evaluate_energy(r);
                } else {
                    auto const pot = derived().sphere_outward_potential(system, i);
                    sum += pot.evaluate_energy(r);
                }
            }

            return sum;
        }

        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::sphere const sphere = derived().sphere_surface(system);
            md::scalar const radius2 = sphere.radius * sphere.radius;

            md::array_view<md::point const> positions = system.view_positions();

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = positions[i] - sphere.center;
                md::scalar const r2 = r.squared_norm();
                md::scalar const r1 = std::sqrt(r2);

                md::vector force;

                if (r2 < radius2) {
                    auto const pot = derived().sphere_inward_potential(system, i);
                    force = pot.evaluate_force(r);
                } else {
                    auto const pot = derived().sphere_outward_potential(system, i);
                    force = pot.evaluate_force(r);
                }

                md::scalar const scale = sphere.radius / r1;
                md::vector const aniso = scale * (force.project(r) - force);

                forces[i] = force + aniso;
            }
        }

        // Default implementations

        md::constant_potential sphere_inward_potential(md::system const&, md::index) const
        {
            return md::constant_potential{0};
        }

        md::constant_potential sphere_outward_potential(md::system const&, md::index) const
        {
            return md::constant_potential{0};
        }

    private:
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }
    };
}

#endif
