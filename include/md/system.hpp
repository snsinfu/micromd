// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_HPP
#define MD_SYSTEM_HPP

#include <algorithm>
#include <memory>

#include "forcefield.hpp"
#include "typedef.hpp"
#include "detail/attribute_table.hpp"
#include "detail/sum_forcefield.hpp"
#include "vendor/array_view.hpp"
#include "vendor/point.hpp"


namespace md
{
    // mass_attribute is an attribute key for particle mass. The default value
    // is 1.
    inline md::scalar mass_attribute(struct tag_mass_attribute*)
    {
        return 1;
    }

    // position_attribute is an attribute key for particle position. The default
    // value is the origin.
    inline md::point position_attribute(struct tag_position_attribute*)
    {
        return {};
    }

    // velocity_attribute is an attribute key for particle velocity. The default
    // value is the zero vector.
    inline md::vector velocity_attribute(struct tag_velocity_attribute*)
    {
        return {};
    }

    // system
    class system
    {
    public:
        system()
        {
            require(md::mass_attribute);
            require(md::position_attribute);
            require(md::velocity_attribute);
        }

        // add_particle adds a particle to the system.
        void add_particle()
        {
            attributes_.resize(attributes_.size() + 1);
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
