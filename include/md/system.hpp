// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_HPP
#define MD_SYSTEM_HPP

#include <algorithm>
#include <memory>

#include "basic_types.hpp"
#include "forcefield.hpp"

#include "system/attribute.hpp"
#include "system/detail/attribute_table.hpp"
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

    //
    struct basic_particle_data
    {
        md::scalar mass = md::default_value(md::mass_attribute);
        md::scalar mobility = md::default_value(md::mobility_attribute);
        md::point position = md::default_value(md::position_attribute);
        md::vector velocity = md::default_value(md::velocity_attribute);
    };

    // system
    class system
    {
    public:
        system()
        {
            require(md::mass_attribute);
            require(md::mobility_attribute);
            require(md::position_attribute);
            require(md::velocity_attribute);
        }

        // add_particle adds a particle to the system.
        void add_particle(md::basic_particle_data data = {})
        {
            md::index const idx = attributes_.size();

            attributes_.resize(idx + 1);

            view(md::mass_attribute)[idx] = data.mass;
            view(md::mobility_attribute)[idx] = data.mobility;
            view(md::position_attribute)[idx] = data.position;
            view(md::velocity_attribute)[idx] = data.velocity;
        }

        // particle_count returns the number of particles in the system.
        md::index particle_count() const
        {
            return attributes_.size();
        }

        // require creates a particle attribute if it does not exist.
        template<typename T, typename Tag>
        void require(md::attribute_key<T, Tag> key)
        {
            attributes_.require(key);
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
        void add_forcefield(std::shared_ptr<md::forcefield> ff)
        {
            forcefield_.add(ff);
        }

        // compute_potential_energy returns the total potential energy of the
        // system.
        md::scalar compute_potential_energy()
        {
            return forcefield_.compute_energy(*this);
        }

        // compute_force assigns total force acting on each particle to given
        // forces array. The length of the array must be the same as the number
        // of particles in the system.
        void compute_force(md::array_view<md::vector> forces)
        {
            assert(forces.size() == particle_count());

            // Zero clear
            std::fill(forces.begin(), forces.end(), md::vector{});

            forcefield_.compute_force(*this, forces);
        }

    private:
        detail::attribute_table attributes_;
        detail::sum_forcefield forcefield_;
    };
}

#endif
