// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_HPP
#define MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_HPP

#include <algorithm>
#include <cmath>
#include <iterator>
#include <utility>
#include <vector>

#include "../../basic_types.hpp"

#include "neighbor_searcher.hpp"


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
            // among heuristics I have tried.
            hash.modulus = md::linear_hash::hash_t(points.size() * 2 / 11);
            hash.modulus |= 1;

            return hash;
        }
    }

    // neighbor_list is a data structure for keeping track of neighbor pairs in
    // a slowly moving particle system.
    class neighbor_list
    {
    public:
        using iterator = std::vector<std::pair<md::index, md::index>>::const_iterator;

        void update(md::array_view<md::point const> points, md::scalar dcut)
        {
            if (!check_consistency(points, dcut)) {
                rebuild(points, dcut);
            }
        }

        iterator begin() const
        {
            return pairs_.begin();
        }

        iterator end() const
        {
            return pairs_.end();
        }

    private:
        bool check_consistency(md::array_view<md::point const> points, md::scalar dcut) const
        {
            if (points.size() != cached_points_.size()) {
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
                if (md::squared_distance(points[i], cached_points_[i]) > threshold * threshold) {
                    return false;
                }
            }

            return true;
        }

        void rebuild(md::array_view<md::point const> points, md::scalar dcut)
        {
            // Rule-of-thumb parameters.
            constexpr md::scalar skin_factor = 1.2;

            verlet_radius_ = dcut * skin_factor;
            cached_points_.assign(points.begin(), points.end());
            pairs_.clear();

            md::linear_hash hash = detail::determine_hash(points, dcut);
            md::neighbor_searcher searcher{verlet_radius_, hash};
            searcher.set_points(cached_points_);
            searcher.search(std::back_inserter(pairs_));
        }

    private:
        md::scalar verlet_radius_ = 0;
        std::vector<md::point> cached_points_;
        std::vector<std::pair<md::index, md::index>> pairs_;
    };
}

#endif
