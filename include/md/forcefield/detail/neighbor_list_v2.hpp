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
#include <utility>
#include <vector>

#include "../../basic_types.hpp"
#include "../../misc/neighbor_searcher_v2.hpp"
#include "neighbor_list_heuristics.hpp"


namespace md
{
    // neighbor_list_v2 is a data structure for efficiently keeping track of
    // neighbor pairs in a slowly moving particle system.
    template<typename Box>
    class neighbor_list_v2
    {
        using iterator = std::vector<std::pair<md::index, md::index>>::const_iterator;

    public:
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
            md::scalar const threshold = (verlet_radius_ - dcut) / 2;

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
            // Let v be the verlet factor. The cost of list construction scales
            // with v^3, while the benefit of list reuse scales with (v-1). So
            // the overall cost scales with v^3/(v-1). The minimum cost is then
            // achieved when v = 1.5.
            static constexpr md::scalar verlet_factor = 1.5;

            verlet_radius_ = verlet_factor * dcut;
            prev_points_.assign(points.begin(), points.end());
            pairs_.clear();

            detail::set_box_hints(box, points);
            md::neighbor_searcher_v2<Box> searcher{box, verlet_radius_};
            searcher.set_points(prev_points_);
            searcher.search(std::back_inserter(pairs_));
        }

    private:
        md::scalar verlet_radius_ = 0;
        std::vector<md::point> prev_points_;
        std::vector<std::pair<md::index, md::index>> pairs_;
    };
}

#endif
