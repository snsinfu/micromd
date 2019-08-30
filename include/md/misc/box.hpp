// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_MISC_BOX_HPP
#define MD_MISC_BOX_HPP

// This module provides box types that represent boundary conditions.

#include "../basic_types.hpp"

#include "nsearch_detail/math.hpp"


namespace md
{
    // open_box represents an open simulation system. No periodic boundaries.
    struct open_box
    {
        // Approximate number of particles in the system. This parameter affects
        // the speed of neighbor search.
        md::index particle_count = 1000;

        // Returns the shortest displacement vector pointing from p2 to p1.
        md::vector shortest_displacement(md::point p1, md::point p2) const
        {
            return p1 - p2;
        }
    };

    // periodic_box represents a system periodic along the coordinate axes.
    struct periodic_box
    {
        // Period along the x axis.
        md::scalar x_period = 1;

        // Period along the y axis.
        md::scalar y_period = 1;

        // Period along the z axis.
        md::scalar z_period = 1;

        // Returns the shortest displacement vector pointing from p2 to p1
        // taking into account the periodicity.
        md::vector shortest_displacement(md::point p1, md::point p2) const
        {
            auto const dx = nsearch_detail::round_mod(p1.x - p2.x, x_period);
            auto const dy = nsearch_detail::round_mod(p1.y - p2.y, y_period);
            auto const dz = nsearch_detail::round_mod(p1.z - p2.z, z_period);
            return {dx, dy, dz};
        }
    };

    // xy_periodic_box represents a system periodic along the x and y axes and
    // open in the z direction.
    struct xy_periodic_box
    {
        // Period along the x axis.
        md::scalar x_period = 1;

        // Period along the y axis.
        md::scalar y_period = 1;

        // Approximate span of the point cloud in the z direction. This
        // parameter affects the speed of neighbor search.
        md::scalar z_span = 1;

        // Approximate number of particles in the system. This parameter affects
        // the speed of neighbor search.
        md::index particle_count = 1000;

        // Returns the shortest displacement vector pointing from p2 to p1
        // taking into account the periodicity.
        md::vector shortest_displacement(md::point p1, md::point p2) const
        {
            auto const dx = nsearch_detail::round_mod(p1.x - p2.x, x_period);
            auto const dy = nsearch_detail::round_mod(p1.y - p2.y, y_period);
            auto const dz = p1.z - p2.z;
            return {dx, dy, dz};
        }
    };
}

#endif
