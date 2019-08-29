// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_MISC_NEIGHBOR_SEARCHER_V2_HPP
#define MD_MISC_NEIGHBOR_SEARCHER_V2_HPP

// This module implements a time-efficient neighbor search algorithm for point
// cloud in open and periodic systems.

#include <algorithm>
#include <utility>

#include "../basic_types.hpp"

#include "nsearch_detail/math.hpp"
#include "nsearch_detail/search_grid.hpp"


namespace md
{
    template<typename Box>
    class neighbor_searcher_v2
    {
        using grid_type = nsearch_detail::search_grid<Box>;

    public:
        // Constructor initializes search strategy using given information.
        // Points are not set.
        neighbor_searcher_v2(Box box, md::scalar dcut)
            : box_{box}, grid_{box, dcut}, dcut_{dcut}
        {
        }

        // Sets points to search.
        void set_points(md::array_view<md::point const> points)
        {
            for (auto& bucket : grid_.buckets) {
                bucket.members.clear();
            }

            for (md::index index = 0; index < points.size(); index++) {
                auto const point = points[index];
                auto const bucket_index = grid_.locate_bucket(point);
                auto& bucket = grid_.buckets[bucket_index];
                bucket.members.push_back({ index, point });
            }
        }

        // Searches neighboring points. Outputs pairs of indices of neighboring
        // points to given output iterator. Each index pair (i, j) satisfies
        // `i < j`. No duplicates are reported.
        template<typename OutputIterator>
        void search(OutputIterator out) const
        {
            for (auto const& bucket : grid_.buckets) {
                for (auto const neighbor_index : bucket.directed_neighbors) {
                    search_among(bucket, grid_.buckets[neighbor_index], out);
                }
            }
        }

        // Searches neighboring points of given one.
        template<typename OutputIterator>
        void query(md::point point, OutputIterator out) const
        {
            auto const bucket_index = grid_.locate_bucket(point);
            auto const& bucket = grid_.buckets[bucket_index];
            query_around(bucket, point, out);
        }

    private:
        // Searches a pair of buckets for neighboring pairs of points.
        template<typename OutputIterator>
        inline void search_among(
            nsearch_detail::spatial_bucket const& bucket_a,
            nsearch_detail::spatial_bucket const& bucket_b,
            OutputIterator& out
        ) const
        {
            auto const dcut2 = dcut_ * dcut_;

            for (auto const member_j : bucket_b.members) {
                for (auto const member_i : bucket_a.members) {
                    if (member_i.index == member_j.index) {
                        // Avoid double counting when bucket_a == bucket_b.
                        break;
                    }

                    if (squared_distance(member_i.point, member_j.point) > dcut2) {
                        continue;
                    }

                    *out++ = std::make_pair(
                        std::min(member_i.index, member_j.index),
                        std::max(member_i.index, member_j.index)
                    );
                }
            }
        }

        // Query neighboring points around the given bucket.
        template<typename OutputIterator>
        inline void query_around(
            nsearch_detail::spatial_bucket const& bucket,
            md::point point,
            OutputIterator& out
        ) const
        {
            auto const dcut2 = dcut_ * dcut_;

            for (auto const neighbor_index : bucket.complete_neighbors) {
                auto const& neighbor = grid_.buckets[neighbor_index];
                for (auto const member : neighbor.members) {
                    if (squared_distance(member.point, point) < dcut2) {
                        *out++ = member.index;
                    }
                }
            }
        }

        inline md::scalar squared_distance(md::point p1, md::point p2) const
        {
            return box_.shortest_displacement(p1, p2).squared_norm();
        }

    private:
        Box box_;
        grid_type grid_;
        md::scalar dcut_;
    };
}

#endif
