// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_INTER_SUBSYSTEM_NEIGHBOR_LIST_HPP
#define MD_FORCEFIELD_DETAIL_INTER_SUBSYSTEM_NEIGHBOR_LIST_HPP

// This module provides inter_subsystem_neighbor_list: A data structure for
// tracking neighbor pairs in two subsystems.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <utility>
#include <vector>

#include "../../basic_types.hpp"
#include "../../misc/linear_hash.hpp"
#include "../../misc/neighbor_searcher.hpp"

#include "neighbor_list_heuristics.hpp"

#include <iostream>

namespace md
{
    // inter_subsystem_neighbor_list is a data structure for efficiently keeping
    // track of neighbor pairs in a subsystem of a slowly moving particle
    // system.
    class inter_subsystem_neighbor_list
    {
        using iterator = std::vector<std::pair<md::index, md::index>>::const_iterator;

    public:
        // set_subsystems sets the subsystems to operate on.
        void set_subsystems(
            md::array_view<md::index const> indices1,
            md::array_view<md::index const> indices2
        )
        {
            key_indices_.assign(indices1.begin(), indices1.end());
            query_indices_.assign(indices2.begin(), indices2.end());

            // Invalidate cache.
            key_points_.clear();
            query_points_.clear();
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

            auto check_within = [&](auto const& sub_indices, auto const& sub_points) {
                if (sub_indices.size() != sub_points.size()) {
                    return false;
                }

                for (std::size_t i = 0; i < sub_indices.size(); i++) {
                    md::index const idx = sub_indices[i];
                    md::point const pt = sub_points[i];
                    if (md::squared_distance(points[idx], pt) > threshold * threshold) {
                        return false;
                    }
                }

                return true;
            };

            if (!check_within(key_indices_, key_points_)) {
                return false;
            }

            if (!check_within(query_indices_, query_points_)) {
                return false;
            }

            return true;
        }

        // rebuild completely rebuilds the neighbor list.
        void rebuild(md::array_view<md::point const> points, md::scalar dcut)
        {
            // Load subsystem points.
            key_points_.clear();
            for (md::index idx : key_indices_) {
                key_points_.push_back(points[idx]);
            }

            query_points_.clear();
            for (md::index idx : query_indices_) {
                query_points_.push_back(points[idx]);
            }

            assert(key_points_.size() == key_indices_.size());
            assert(query_points_.size() == query_indices_.size());

            // Hash first subsystem and query second subsystem for neighbors.
            verlet_radius_ = detail::determine_verlet_radius(dcut);
            md::linear_hash hash = detail::determine_hash(points, verlet_radius_);
            md::neighbor_searcher searcher{verlet_radius_, hash};
            searcher.set_points(key_points_);

            pairs_.clear();
            std::vector<md::index> hits;

            for (std::size_t i = 0; i < query_points_.size(); i++) {
                md::point const query_pt = query_points_[i];
                md::index const query_idx = query_indices_[i];

                hits.clear();
                searcher.query(query_pt, std::back_inserter(hits));

                for (md::index const j : hits) {
                    md::index const key_idx = key_indices_[j];
                    pairs_.push_back(
                        std::make_pair(
                            std::min(key_idx, query_idx),
                            std::max(key_idx, query_idx)
                        )
                    );
                }
            }
        }

    private:
        md::scalar verlet_radius_ = 0;
        std::vector<md::index> key_indices_;
        std::vector<md::index> query_indices_;
        std::vector<md::point> key_points_;
        std::vector<md::point> query_points_;
        std::vector<std::pair<md::index, md::index>> pairs_;
    };
}

#endif
