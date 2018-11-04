// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_HPP
#define MD_FORCEFIELD_DETAIL_NEIGHBOR_LIST_HPP

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include "../../basic_types.hpp"

#include "neighbor_searcher.hpp"


namespace md
{
    namespace detail
    {
        // compute_variance returns the population variance of each coordinate
        // of given points as a vector. Returns a zero vector for empty input.
        inline md::vector compute_variance(md::array_view<md::point const> points)
        {
            if (points.empty()) {
                return {};
            }

            md::scalar sum_x = 0;
            md::scalar sum_y = 0;
            md::scalar sum_z = 0;

            md::scalar sum_x2 = 0;
            md::scalar sum_y2 = 0;
            md::scalar sum_z2 = 0;

            for (md::point const& pt : points) {
                sum_x += pt.x;
                sum_y += pt.y;
                sum_z += pt.z;

                sum_x2 += pt.x * pt.x;
                sum_y2 += pt.y * pt.y;
                sum_z2 += pt.z * pt.z;
            }

            md::scalar const n = md::scalar(points.size());

            return md::vector {
                sum_x2 / n - (sum_x / n) * (sum_x / n),
                sum_y2 / n - (sum_y / n) * (sum_y / n),
                sum_z2 / n - (sum_z / n) * (sum_z / n)
            };
        }

        // determine_hash returns a linear_hash object that is heuristically
        // parameterized to make neighbor_searcher perform good on given points.
        inline md::linear_hash determine_hash(md::array_view<md::point const> points)
        {
            // Heuristic: Squash lowest-variance axis.
            md::vector const var = detail::compute_variance(points);
            md::scalar const min = std::min({var.x, var.y, var.z});

            md::linear_hash hash;

            if (var.x == min) {
                hash.x_coeff = 0;
            } else if (var.y == min) {
                hash.y_coeff = 0;
            } else {
                hash.z_coeff = 0;
            }

            // Heuristic: 10-30 give good bucket distribution for dense system.
            constexpr md::index bucket_occupancy = 16;

            hash.modulus = md::linear_hash::hash_t(points.size() / bucket_occupancy);
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

            md::linear_hash hash = detail::determine_hash(points);
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
