// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_INTRA_SUBSYSTEM_NEIGHBOR_LIST_HPP
#define MD_FORCEFIELD_DETAIL_INTRA_SUBSYSTEM_NEIGHBOR_LIST_HPP

// This module provides intra_subsystem_neighbor_list: A data structure for
// tracking neighbor pairs in a subsystem.

#include <algorithm>
#include <cassert>
#include <iterator>
#include <utility>
#include <vector>

#include "../../basic_types.hpp"
#include "../../misc/linear_hash.hpp"
#include "../../misc/neighbor_searcher.hpp"

#include "neighbor_list_heuristics.hpp"


namespace md
{
    // intra_subsystem_neighbor_list is a data structure for efficiently keeping
    // track of neighbor pairs in a subsystem of a slowly moving particle
    // system.
    class intra_subsystem_neighbor_list
    {
        using iterator = std::vector<std::pair<md::index, md::index>>::const_iterator;

    public:
        // set_subsystem sets the subsystem to operate on.
        void set_subsystem(md::array_view<md::index const> indices)
        {
            subsystem_.clear();
            for (md::index idx : indices) {
                point_index pt;
                pt.index = idx;
                subsystem_.push_back(pt);
            }
        }

        // update rebuilds the neighbor list if necessary.
        void update(md::array_view<md::point const> points, md::scalar dcut)
        {
            if (!check_consistency(points, dcut)) {
                rebuild(points, dcut);
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
        // check_consistency checks if the cached neighbor list is still usable
        // with given points and cutoff distance.
        bool check_consistency(md::array_view<md::point const> points, md::scalar dcut) const
        {
            // False negatives (unlisted point pairs that fall actually within
            // dcut) won't arise if the displacement from previous rebuild is
            // less than or equal to this threshold.
            md::scalar const threshold = (verlet_radius_ - dcut) / 2;

            if (threshold <= 0) {
                return false;
            }

            for (point_index const& pt : subsystem_) {
                assert(pt.index < points.size());
                if (md::squared_distance(points[pt.index], pt.point) > threshold * threshold) {
                    return false;
                }
            }

            return true;
        }

        // rebuild completely rebuilds the neighbor list.
        void rebuild(md::array_view<md::point const> points, md::scalar dcut)
        {
            for (point_index& pt : subsystem_) {
                assert(pt.index < points.size());
                pt.point = points[pt.index];
            }

            // FIXME: Dumb copying.
            std::vector<md::point> subsystem_points;
            subsystem_points.reserve(subsystem_.size());
            for (point_index const& pt : subsystem_) {
                subsystem_points.push_back(pt.point);
            }

            pairs_.clear();
            verlet_radius_ = detail::determine_verlet_radius(dcut);
            md::linear_hash hash = detail::determine_hash(points, verlet_radius_);
            md::neighbor_searcher searcher{verlet_radius_, hash};
            searcher.set_points(subsystem_points);
            searcher.search(std::back_inserter(pairs_));

            // Translate back the indices.
            for (auto& pair : pairs_) {
                pair.first = subsystem_[pair.first].index;
                pair.second = subsystem_[pair.second].index;
            }
        }

    private:
        struct point_index
        {
            md::point point;
            md::index index;
        };

        md::scalar verlet_radius_ = 0;
        std::vector<point_index> subsystem_;
        std::vector<std::pair<md::index, md::index>> pairs_;
    };
}

#endif
