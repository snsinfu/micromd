// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_DETAIL_NEIGHBOR_SEARCHER_HPP
#define MD_DETAIL_NEIGHBOR_SEARCHER_HPP

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include "linear_hash.hpp"
#include "../typedef.hpp"
#include "../vendor/array_view.hpp"
#include "../vendor/point.hpp"


namespace md
{
    // neighbor_searcher
    class neighbor_searcher
    {
    public:
        // Constructor takes the cut-off distance and a hash function.
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

            for (md::index center = 0; center < hash.modulus; center++) {
                std::vector<md::index>& neighbors = buckets_[center].neighbors;

                for (hash_t const delta : hash_deltas) {
                    hash_t const neighbor = (center + delta) % hash.modulus;

                    // Leverage symmetry.
                    if (neighbor >= center) {
                        neighbors.push_back(neighbor);
                    }
                }
            }
        }

        // set_points fills hash table.
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

        // search
        template<typename OutputIterator>
        void search(OutputIterator out)
        {
            for (hash_bucket const& center : buckets_) {
                for (md::index const neighbor : center.neighbors) {
                    search_among(center, buckets_[neighbor], out);
                }
            }
        }

    private:
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

        template<typename OutputIterator>
        inline void search_among(
            hash_bucket const& bucket_a,
            hash_bucket const& bucket_b,
            OutputIterator& out
        )
        {
            md::scalar const dcut2 = dcut_ * dcut_;

            for (hash_bucket::member member_j : bucket_b.members) {
                for (hash_bucket::member member_i : bucket_a.members) {
                    if (member_i.index == member_j.index) {
                        // Avoid double counting in intra-bucket search.
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

            // Offset coordinates to avoid negative value.
            constexpr md::scalar offset = 1L << 30;

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
