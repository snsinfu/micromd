// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_V2_HPP
#define MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_V2_HPP

// This internal module provides neighbor_list_v2: A data structure for tracking
// neighbor pairs in a system. Used to implement neighbor_pair_forcefield_v2.

#include <algorithm>
#include <cmath>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

#include "../../basic_types.hpp"
#include "../../misc/box.hpp"
#include "../../misc/neighbor_searcher_v2.hpp"
#include "neighbor_list_heuristics.hpp"


namespace md
{
    namespace detail
    {
        inline bool approx(md::scalar x, md::scalar y)
        {
            // FIXME: Magic number
            constexpr md::scalar epsilon = 1e-6;
            return std::fabs(x - y) < epsilon * std::fabs(y);
        }

        // FIXME: Box abstraction leaks here
        inline bool approx(md::open_box box1, md::open_box box2)
        {
            return box1.particle_count == box2.particle_count;
        }

        inline bool approx(md::periodic_box box1, md::periodic_box box2)
        {
            return approx(box1.x_period, box2.x_period)
                && approx(box1.y_period, box2.y_period)
                && approx(box1.z_period, box2.z_period);
        }

        inline bool approx(md::xy_periodic_box box1, md::xy_periodic_box box2)
        {
            return approx(box1.x_period, box2.x_period)
                && approx(box1.y_period, box2.y_period)
                && approx(box1.z_span, box2.z_span)
                && box1.particle_count == box2.particle_count;
        }
    }

    // neighbor_list_v2 is a data structure for efficiently keeping track of
    // neighbor pairs in a slowly moving particle system.
    template<typename Box>
    class neighbor_list_v2
    {
        using iterator = std::vector<std::pair<md::index, md::index>>::const_iterator;

    public:
        neighbor_list_v2()
            : searcher_{prev_box_, prev_verlet_radius_} // FIXME: Poor default
        {
        }

        // Rebuilds the neighbor list if necessary.
        void update(md::array_view<md::point const> points, md::scalar dcut, Box box)
        {
            if (!check_consistency(points, dcut, box)) {
                rebuild(points, dcut, box);
            }
        }

        // Range interface.
        iterator begin() const
        {
            return pairs_.begin();
        }

        iterator end() const
        {
            return pairs_.end();
        }

    private:
        // Checks if the previously created neighbor list is still usable with
        // the given configuration.
        bool check_consistency(
            md::array_view<md::point const> points, md::scalar dcut, Box box
        ) const
        {
            if (points.size() != prev_points_.size()) {
                return false;
            }

            // False negatives (unlisted point pairs that fall actually within
            // dcut) won't arise if the displacement from previous rebuild is
            // less than or equal to this threshold.
            md::scalar const threshold = (prev_verlet_radius_ - dcut) / 2;

            if (threshold <= 0) {
                return false;
            }

            for (md::index i = 0; i < points.size(); i++) {
                md::vector const disp = box.shortest_displacement(points[i], prev_points_[i]);
                if (disp.squared_norm() > threshold * threshold) {
                    return false;
                }
            }
            return true;
        }

        // Rebuilds the neighbor list.
        void rebuild(
            md::array_view<md::point const> points, md::scalar dcut, Box box
        )
        {
            detail::set_box_hints(box, points);

            // Let v be the verlet factor. The cost of list construction scales
            // with v^3, while the benefit of list reuse scales with (v-1). So
            // the overall cost scales with v^3/(v-1). The minimum cost is then
            // achieved when v = 1.5.
            static constexpr md::scalar verlet_factor = 1.5;
            md::scalar const verlet_radius = verlet_factor * dcut;

            // Neighbor searcher is expensive to construct. Reuse previous one
            // if possible.
            bool const box_changed = !detail::approx(box, prev_box_);
            bool const verlet_changed = !detail::approx(verlet_radius, prev_verlet_radius_);
            if (box_changed || verlet_changed) {
                searcher_ = md::neighbor_searcher_v2<Box>{box, verlet_radius};
            }

            prev_points_.assign(points.begin(), points.end());
            pairs_.clear();
            searcher_.set_points(prev_points_);
            searcher_.search(std::back_inserter(pairs_));
        }

    private:
        Box prev_box_;
        md::scalar prev_verlet_radius_ = 1;
        md::neighbor_searcher_v2<Box> searcher_;
        std::vector<md::point> prev_points_;
        std::vector<std::pair<md::index, md::index>> pairs_;
    };
}

#endif
