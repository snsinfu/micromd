// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_NEIGHBOR_SEARCHER_HPP
#define MD_FORCEFIELD_DETAIL_NEIGHBOR_SEARCHER_HPP

// This module implements a time-efficient neighbor search algorithm in an open
// three-dimensional space.

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include "../../basic_types.hpp"

#include "linear_hash.hpp"


namespace md
{
    // neighbor_searcher implements a spatial hashing algorithm for efficiently
    // searching given points for neighboring pairs.
    class neighbor_searcher
    {
    public:
        // Constructor takes the cut-off distance and a hash function to use.
        //
        // The modulus of the hash determines the number of buckets where points
        // are assigned. It should not be .
        neighbor_searcher(md::scalar dcut, md::linear_hash hash)
            : dcut_{dcut}, hash_{hash}, buckets_(hash.modulus)
        {
            // Pre-compute neighborhood of each bucket for faster query.
            using hash_t = linear_hash::hash_t;

            hash_t const coord_deltas[] = {
                hash.modulus - 1,
                hash.modulus,
                hash.modulus + 1
            };
            std::set<hash_t> hash_deltas;

            for (hash_t const dx : coord_deltas) {
                for (hash_t const dy : coord_deltas) {
                    for (hash_t const dz : coord_deltas) {
                        hash_deltas.insert(hash_(dx, dy, dz));
                    }
                }
            }

            for (hash_t center = 0; center < hash.modulus; center++) {
                std::vector<md::index>& neighbors = buckets_[center].neighbors;

                for (hash_t const delta : hash_deltas) {
                    hash_t const neighbor = (center + delta) % hash.modulus;

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
        inline linear_hash::hash_t locate_bucket(md::point pt) const
        {
            using hash_t = linear_hash::hash_t;

            // Negative coordinate value causes discontinuous jumps in hash
            // value which breaks our search algorithm. Avoid that by
            // offsetting.
            constexpr md::scalar offset = 1L << 16;

            md::scalar const freq = 1 / dcut_;
            hash_t const x = hash_t(offset + freq * pt.x);
            hash_t const y = hash_t(offset + freq * pt.y);
            hash_t const z = hash_t(offset + freq * pt.z);

            return hash_(x, y, z);
        }

    private:
        md::scalar dcut_;
        md::linear_hash hash_;
        std::vector<hash_bucket> buckets_;
    };
}

#endif
