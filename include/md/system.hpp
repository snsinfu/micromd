// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_HPP
#define MD_SYSTEM_HPP

// This module defines the system class: A context class holding all particle
// parameters and state variables.

#include <algorithm>
#include <memory>
#include <type_traits>

#include "basic_types.hpp"
#include "forcefield.hpp"

#include "system/attribute.hpp"
#include "system/particle.hpp"
#include "system/detail/attribute_table.hpp"
#include "system/detail/iterator_range.hpp"
#include "system/detail/sum_forcefield.hpp"


namespace md
{
    // mass_attribute is an attribute key for particle mass. The default value
    // is 1.
    inline constexpr md::scalar mass_attribute(struct tag_mass_attribute*)
    {
        return 1;
    }

    // mobility_attribute is an attribute key for particle mobility. This is
    // used in brownian dynamics simulations. The default value is 1.
    inline constexpr md::scalar mobility_attribute(struct tag_mobility_attribute*)
    {
        return 1;
    }

    // position_attribute is an attribute key for particle position. The default
    // value is the origin.
    inline constexpr md::point position_attribute(struct tag_position_attribute*)
    {
        return {};
    }

    // velocity_attribute is an attribute key for particle velocity. The default
    // value is the zero vector.
    inline constexpr md::vector velocity_attribute(struct tag_velocity_attribute*)
    {
        return {};
    }

    // basic_particle_data holds basic particle attribute values. It is used to
    // pass these data to system::add_particle function.
    struct basic_particle_data
    {
        md::scalar mass = md::default_value(md::mass_attribute);
        md::scalar mobility = md::default_value(md::mobility_attribute);
        md::point position = md::default_value(md::position_attribute);
        md::vector velocity = md::default_value(md::velocity_attribute);
    };

    // system is a context class.
    class system
    {
    public:
        system()
        {
            add_attribute(md::mass_attribute);
            add_attribute(md::mobility_attribute);
            add_attribute(md::position_attribute);
            add_attribute(md::velocity_attribute);
        }

        // add_particle adds a particle to the system.
        md::particle_ref add_particle(md::basic_particle_data data = {})
        {
            md::index const idx = attributes_.size();

            attributes_.resize(idx + 1);

            view_masses()[idx] = data.mass;
            view_mobilities()[idx] = data.mobility;
            view_positions()[idx] = data.position;
            view_velocities()[idx] = data.velocity;

            return md::particle_ref(*this, idx);
        }

        // particles returns a range for iterating over particles in the system.
        detail::iterator_range<md::particle_iterator> particles()
        {
            return detail::iterator_range<md::particle_iterator>{
                md::particle_iterator{*this, 0},
                md::particle_iterator{*this, particle_count()},
            };
        }

        // particle_count returns the number of particles in the system.
        md::index particle_count() const
        {
            return attributes_.size();
        }

        // add_attribute creates a particle attribute if it does not exist.
        template<typename T, typename Tag>
        void add_attribute(md::attribute_key<T, Tag> key)
        {
            attributes_.require(key);
        }

        // `require` is an alias of `add_attribute`. It will be deprecated in
        // the future.
        template<typename T, typename Tag>
        void require(md::attribute_key<T, Tag> key)
        {
            add_attribute(key);
        }

        // view returns a view into the array of particle attribute values.
        template<typename T, typename Tag>
        md::array_view<T> view(md::attribute_key<T, Tag> key)
        {
            return attributes_.view(key);
        }

        template<typename T, typename Tag>
        md::array_view<T const> view(md::attribute_key<T, Tag> key) const
        {
            return attributes_.view(key);
        }

        // view_masses returns a view of built-in mass attributes.
        md::array_view<md::scalar> view_masses() noexcept
        {
            return attributes_.view(md::mass_attribute);
        }

        md::array_view<md::scalar const> view_masses() const noexcept
        {
            return attributes_.view(md::mass_attribute);
        }

        // view_mobilities returns a view of built-in mobility attributes.
        md::array_view<md::scalar> view_mobilities() noexcept
        {
            return attributes_.view(md::mobility_attribute);
        }

        md::array_view<md::scalar const> view_mobilities() const noexcept
        {
            return attributes_.view(md::mobility_attribute);
        }

        // view_positions returns a view of built-in position attributes.
        md::array_view<md::point> view_positions() noexcept
        {
            return attributes_.view(md::position_attribute);
        }

        md::array_view<md::point const> view_positions() const noexcept
        {
            return attributes_.view(md::position_attribute);
        }

        // view_velocities returns a view of built-in velocity attributes.
        md::array_view<md::vector> view_velocities() noexcept
        {
            return attributes_.view(md::velocity_attribute);
        }

        md::array_view<md::vector const> view_velocities() const noexcept
        {
            return attributes_.view(md::velocity_attribute);
        }

        // add_forcefield adds a forcefield to the system.
        //
        // With this overload the added forcefield is shared between the caller.
        // So the added forcefield can be modified from outside.
        void add_forcefield(std::shared_ptr<md::forcefield> ff)
        {
            forcefield_.add(ff);
        }

        // add_forcefield adds a copy of passed forcefield to the system.
        //
        // This overload returns a new shared_ptr holding the copy of the passed
        // forcefield. The forcefield can be modified via the returned pointer.
        template<
            typename FF,
            typename = typename std::enable_if<std::is_base_of<md::forcefield, FF>::value>::type
        >
        std::shared_ptr<FF> add_forcefield(FF const& ff)
        {
            std::shared_ptr<FF> ffptr = std::make_shared<FF>(ff);
            add_forcefield(ffptr);
            return ffptr;
        }

        // compute_kinetic_energy returns the total kinetic energy of the
        // system. It uses mass and velocity particle attributes.
        md::scalar compute_kinetic_energy() const
        {
            md::array_view<md::scalar const> masses = view_masses();
            md::array_view<md::vector const> velocities = view_velocities();

            md::scalar sum = 0;

            for (md::index i = 0; i < particle_count(); i++) {
                sum += 0.5 * masses[i] * velocities[i].squared_norm();
            }

            return sum;
        }

        // compute_potential_energy returns the total potential energy of the
        // system. This function is not marked const because forcefield may
        // mutate its state when computing energy.
        md::scalar compute_potential_energy()
        {
            return forcefield_.compute_energy(*this);
        }

        // compute_energy returns the total mechanical energy of the system.
        // Mechanical energy is the sum of kinetic energy and potential energy.
        md::scalar compute_energy()
        {
            return compute_kinetic_energy() + compute_potential_energy();
        }

        // compute_force assigns total force acting on each particle to given
        // forces array. The length of the array must be the same as the number
        // of particles in the system.
        void compute_force(md::array_view<md::vector> forces)
        {
            assert(forces.size() == particle_count());

            md::vector const zero {};
            std::fill(forces.begin(), forces.end(), zero);

            forcefield_.compute_force(*this, forces);
        }

    private:
        detail::attribute_table attributes_;
        detail::sum_forcefield forcefield_;
    };

    inline md::particle_ref::particle_ref(md::system& system, md::index idx)
        : system{system}
        , index{idx}
        , mass{system.view_masses()[idx]}
        , mobility{system.view_mobilities()[idx]}
        , position{system.view_positions()[idx]}
        , velocity{system.view_velocities()[idx]}
    {
    }

    template<typename T, typename Tag>
    T& md::particle_ref::view(md::attribute_key<T, Tag> key)
    {
        return system.view(key)[index];
    }
}

#endif
