// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_MISC_NEIGHBOR_SEARCHER_HPP
#define MD_MISC_NEIGHBOR_SEARCHER_HPP

// This module implements a time-efficient neighbor search algorithm in an open
// three-dimensional space.

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include "../basic_types.hpp"

#include "linear_hash.hpp"


namespace md
{
    // neighbor_searcher implements a spatial hashing algorithm for efficiently
    // searching point cloud for neighboring pairs.
    class neighbor_searcher
    {
        using hash_type = md::linear_hash::uint;

    public:
        // Constructor takes the cut-off distance and a hash function to use.
        //
        // The modulus of the hash determines the number of buckets where points
        // are assigned. It should not be .
        neighbor_searcher(md::scalar dcut, md::linear_hash hash)
            : dcut_{dcut}, hash_{hash}, buckets_(hash.modulus)
        {
            // Pre-compute neighborhood of each bucket for faster query.
            hash_type const coord_deltas[] = {
                hash.modulus - 1,
                hash.modulus,
                hash.modulus + 1
            };

            for (hash_type const dx : coord_deltas) {
                for (hash_type const dy : coord_deltas) {
                    for (hash_type const dz : coord_deltas) {
                        hash_deltas_.insert(hash_(dx, dy, dz));
                    }
                }
            }

            for (hash_type center = 0; center < hash.modulus; center++) {
                std::vector<md::index>& neighbors = buckets_[center].neighbors;

                for (hash_type const delta : hash_deltas_) {
                    hash_type const neighbor = (center + delta) % hash.modulus;

                    // Leverage symmetry to reduce search space.
                    if (neighbor >= center) {
                        neighbors.push_back(neighbor);
                    }
                }
            }
        }

        // set_points assigns given points to hash table for subsequent queries.
        void set_points(md::array_view<md::point const> points)
        {
            for (hash_bucket& bucket : buckets_) {
                bucket.members.clear();
            }

            for (md::index idx = 0; idx < points.size(); idx++) {
                md::point const pt = points[idx];
                hash_bucket& bucket = buckets_[locate_bucket(pt)];
                bucket.members.push_back({idx, pt});
            }
        }

        // search outputs pairs of indices of neighboring points to given output
        // iterator. Each index pair (i,j) satisfies i < j. No duplicates are
        // reported.
        template<typename OutputIterator>
        void search(OutputIterator out) const
        {
            for (hash_bucket const& center : buckets_) {
                for (md::index const neighbor : center.neighbors) {
                    search_among(center, buckets_[neighbor], out);
                }
            }
        }

        // query outputs the indices of points neighboring to given point.
        template<typename OutputIterator>
        void query(md::point point, OutputIterator out) const
        {
            hash_type const center_index = locate_bucket(point);
            md::scalar const dcut2 = dcut_ * dcut_;

            for (hash_type const delta : hash_deltas_) {
                hash_bucket const& bucket = buckets_[(center_index + delta) % hash_.modulus];

                for (hash_bucket::member member : bucket.members) {
                    if (md::squared_distance(member.point, point) < dcut2) {
                        *out++ = member.index;
                    }
                }
            }

        }

    private:
        // hash_bucket is a bucket in a spatial hash table.
        struct hash_bucket
        {
            struct member
            {
                md::index index;
                md::point point;
            };
            std::vector<member> members;
            std::vector<md::index> neighbors;
        };

        // search_among searches a pair of buckets for all neighboring pairs.
        template<typename OutputIterator>
        inline void search_among(
            hash_bucket const& bucket_a,
            hash_bucket const& bucket_b,
            OutputIterator& out
        ) const
        {
            md::scalar const dcut2 = dcut_ * dcut_;

            for (hash_bucket::member member_j : bucket_b.members) {
                for (hash_bucket::member member_i : bucket_a.members) {
                    if (member_i.index == member_j.index) {
                        // Avoid double counting when bucket_a == bucket_b.
                        break;
                    }

                    if (md::squared_distance(member_i.point, member_j.point) > dcut2) {
                        continue;
                    }

                    *out++ = std::make_pair(
                        std::min(member_i.index, member_j.index),
                        std::max(member_i.index, member_j.index)
                    );
                }
            }
        }

        // locate_bin returns the hash index for a point.
        inline hash_type locate_bucket(md::point pt) const
        {
            // Negative coordinate value causes discontinuous jumps in hash
            // value which breaks our search algorithm. Avoid that by
            // offsetting.
            constexpr md::scalar offset = 1L << 16;

            md::scalar const freq = 1 / dcut_;
            hash_type const x = hash_type(offset + freq * pt.x);
            hash_type const y = hash_type(offset + freq * pt.y);
            hash_type const z = hash_type(offset + freq * pt.z);

            return hash_(x, y, z);
        }

    private:
        md::scalar dcut_;
        md::linear_hash hash_;
        std::set<hash_type> hash_deltas_;
        std::vector<hash_bucket> buckets_;
    };
}

#endif