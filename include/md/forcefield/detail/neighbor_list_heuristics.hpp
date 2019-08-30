// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_HEURISTICS_HPP
#define MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_HEURISTICS_HPP

// This module provides some heuristic helper functions for neighbor_list and
// subsystem_neighbor_list implementation.

#include <cmath>

#include "../../basic_types.hpp"
#include "../../misc/box.hpp"
#include "../../misc/linear_hash.hpp"


namespace md
{
    namespace detail
    {
        // determine_hash returns a linear_hash object that is heuristically
        // parameterized to make neighbor_searcher perform good on given points.
        inline md::linear_hash determine_hash(
            md::array_view<md::point const> points,
            md::scalar dcut
        )
        {
            md::linear_hash hash;

            (void) dcut;

            // Benchmark simulations run fastest with this simple heuristic
            // among those I have tried.
            hash.modulus = md::linear_hash::uint(points.size() * 2 / 11);
            hash.modulus |= 1;

            return hash;
        }

        // determine_verlet_radius determines a good vertlet cutoff radius for
        // neighbor list calculation.
        inline md::scalar determine_verlet_radius(md::scalar dcut)
        {
            // Heuristic: 1.2 gives fairly good performance.
            return 1.2 * dcut;
        }

        // stddev_points computes the standard deviation of points along each
        // axis.
        inline md::vector stddev_points(md::array_view<md::point const> points)
        {
            md::vector mean;
            md::vector mean_sq;

            for (auto const point : points) {
                auto const r = point - points[0];
                mean += r;
                mean_sq += r.hadamard(r);
            }

            if (!points.empty()) {
                mean /= md::scalar(points.size());
                mean_sq /= md::scalar(points.size());
            }

            auto const var = mean_sq - mean.hadamard(mean);
            auto const std_x = std::sqrt(var.x);
            auto const std_y = std::sqrt(var.y);
            auto const std_z = std::sqrt(var.z);
            return {std_x, std_y, std_z};
        }

        // set_box_hints updates hint fields (e.g., z_span) of box using some
        // heuristics on given points.
        template<typename Box>
        void set_box_hints(
            Box& box, md::array_view<md::point const> points
        )
        {
            // No hint in general case.
            (void) box;
            (void) points;
        }

        template<>
        inline void set_box_hints<md::open_box>(
            md::open_box& box, md::array_view<md::point const> points
        )
        {
            box.particle_count = points.size();
        }

        template<>
        inline void set_box_hints<md::xy_periodic_box>(
            md::xy_periodic_box& box, md::array_view<md::point const> points
        )
        {
            box.particle_count = points.size();

            // The span of uniform distribution is about 3.5x the stddev. The
            // real distribution is not necessarily uniform but it works well
            // in practice.
            constexpr md::scalar span_per_stddev = 3.5;
            auto const stddev = stddev_points(points);
            box.z_span = span_per_stddev * stddev.z;
        }
    }
}

#endif
