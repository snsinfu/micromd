// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_MISC_NSEARCH_DETAIL_SEARCH_GRID_HPP
#define MD_MISC_NSEARCH_DETAIL_SEARCH_GRID_HPP

// This module defines spatial_bucket, and search_grid for standard boxes.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "../../basic_types.hpp"
#include "../box.hpp"
#include "../linear_hash.hpp"
#include "math.hpp"


namespace md
{
    namespace nsearch_detail
    {
        using std::uint32_t;
        using std::size_t;

        // Sort and remove duplicates in a container.
        template<typename Container>
        void sort_unique(Container& cont)
        {
            std::sort(cont.begin(), cont.end());
            cont.erase(std::unique(cont.begin(), cont.end()), cont.end());
        }

        // Bucket for spatial hash table.
        struct spatial_bucket
        {
            struct member
            {
                md::index index;
                md::point point;
            };

            // Points in the bucket.
            std::vector<member> members;

            // Indices of the buckets spatially adjoining to this bucket.
            std::vector<md::index> complete_neighbors;

            // Subset of complete_neighbors where bidirectional links in the
            // adjoining graph are reduced to a unidirectional one.
            std::vector<md::index> directed_neighbors;
        };

        // Helper class for generating x/y/z_bins for given box and spacing.
        template<typename Box>
        struct basic_binner;

        // search_grid abstracts away how points should be grouped into spatial
        // neighbors in given box.
        template<typename Box>
        struct search_grid;

        // Generic implementation of search_grid for bin-based construction.
        template<typename Box>
        struct binned_search_grid
        {
            nsearch_detail::bin_layout  x_bins;
            nsearch_detail::bin_layout  y_bins;
            nsearch_detail::bin_layout  z_bins;
            std::vector<spatial_bucket> buckets;

            binned_search_grid(Box box, md::scalar spacing)
            {
                init_bins(box, spacing);
                init_buckets();
            }

            inline size_t locate_bucket(md::point pt) const
            {
                auto const x = nsearch_detail::locate_bin(x_bins, pt.x);
                auto const y = nsearch_detail::locate_bin(y_bins, pt.y);
                auto const z = nsearch_detail::locate_bin(z_bins, pt.z);
                return do_locate_bucket(x, y, z);
            }

        private:
            void init_bins(Box box, md::scalar spacing)
            {
                basic_binner<Box> binner{box, spacing};
                x_bins = binner.x_bins;
                y_bins = binner.y_bins;
                z_bins = binner.z_bins;
            }

            void init_buckets()
            {
                buckets.clear();
                buckets.resize(x_bins.count * y_bins.count * z_bins.count);

                for (uint32_t z = 0; z < z_bins.count; z++) {
                    for (uint32_t y = 0; y < y_bins.count; y++) {
                        for (uint32_t x = 0; x < x_bins.count; x++) {
                            init_bucket_adjacents(x, y, z);
                        }
                    }
                }
            }

            // Initializes adjacent links of a bucket that covers (x,y,z) bin.
            void init_bucket_adjacents(uint32_t x, uint32_t y, uint32_t z)
            {
                auto const bucket_index = do_locate_bucket(x, y, z);
                auto& bucket = buckets[bucket_index];

                uint32_t const dx_values[] = {0, 1, x_bins.count - 1};
                uint32_t const dy_values[] = {0, 1, y_bins.count - 1};
                uint32_t const dz_values[] = {0, 1, z_bins.count - 1};

                for (auto const dz : dz_values) {
                    for (auto const dy : dy_values) {
                        for (auto const dx : dx_values) {
                            auto const adj_x = nsearch_detail::add_mod(x, dx, x_bins.count);
                            auto const adj_y = nsearch_detail::add_mod(y, dy, y_bins.count);
                            auto const adj_z = nsearch_detail::add_mod(z, dz, z_bins.count);
                            auto const adj_index = do_locate_bucket(adj_x, adj_y, adj_z);
                            if (adj_index <= bucket_index) {
                                bucket.directed_neighbors.push_back(adj_index);
                            }
                            bucket.complete_neighbors.push_back(adj_index);
                        }
                    }
                }

                sort_unique(bucket.directed_neighbors);
                sort_unique(bucket.complete_neighbors);
            }

            size_t do_locate_bucket(uint32_t x, uint32_t y, uint32_t z) const
            {
                return x + x_bins.count * (y + y_bins.count * z);
            }
        };

        // Periodic box is easy. Just split each axis into uniform bins.
        template<>
        struct basic_binner<md::periodic_box>
        {
            bin_layout x_bins;
            bin_layout y_bins;
            bin_layout z_bins;

            basic_binner(md::periodic_box box, md::scalar spacing)
            {
                x_bins = nsearch_detail::define_bins(box.x_period, spacing);
                y_bins = nsearch_detail::define_bins(box.y_period, spacing);
                z_bins = nsearch_detail::define_bins(box.z_period, spacing);
            }
        };

        template<>
        struct search_grid<md::periodic_box> : binned_search_grid<md::periodic_box>
        {
            using binned_search_grid::binned_search_grid;
        };

        // xy-periodic system. We try to define the bins along the z direction
        // so that the resulting buckets do not get too sparse.
        template<>
        struct basic_binner<md::xy_periodic_box>
        {
            nsearch_detail::bin_layout x_bins;
            nsearch_detail::bin_layout y_bins;
            nsearch_detail::bin_layout z_bins;

            basic_binner(md::xy_periodic_box box, md::scalar spacing)
            {
                x_bins = nsearch_detail::define_bins(box.x_period, spacing);
                y_bins = nsearch_detail::define_bins(box.y_period, spacing);
                z_bins = estimate_z_bins(box, spacing);
            }

        private:
            nsearch_detail::bin_layout estimate_z_bins(
                md::xy_periodic_box box, md::scalar spacing
            ) const
            {
                // Too low bucket occupancy makes neighbor search inefficient.
                // This is a heuristic minimum.
                constexpr md::scalar min_bin_occupancy = 4.0;

                const auto volume = box.x_period * box.y_period * box.z_span;
                const auto density = md::scalar(box.particle_count) / volume;
                const auto bin_volume = x_bins.step * y_bins.step * spacing;
                const auto bin_occupancy = density * bin_volume;

                const auto target_occupancy = std::max(bin_occupancy, min_bin_occupancy);
                const auto bucket_count = md::scalar(box.particle_count) / target_occupancy;
                const auto z_mod = bucket_count / (x_bins.count * y_bins.count);

                nsearch_detail::bin_layout bins;
                bins.step = spacing;
                bins.count = nsearch_detail::round_uint(z_mod < 1 ? 1 : z_mod);
                return bins;
            }
        };

        template<>
        struct search_grid<md::xy_periodic_box> : binned_search_grid<md::xy_periodic_box>
        {
            using binned_search_grid::binned_search_grid;
        };

        // search_grid implementation for open_box. Open system tends to be
        // sparse, so we use hashing instead of binning to construct a grid.
        template<>
        struct search_grid<md::open_box>
        {
            md::scalar                  spacing;
            md::linear_hash             hash;
            std::vector<spatial_bucket> buckets;

            search_grid(md::open_box box, md::scalar spacing)
                : spacing{spacing}
            {
                init_hash(box);
                init_buckets();
            }

            inline size_t locate_bucket(md::point pt) const
            {
                // Negative coordinate value causes discontinuous jumps in hash value
                // which breaks our search algorithm. Avoid that by offsetting. The
                // offset should be chosen to not overwhelm the coordinate values.
                constexpr md::scalar offset = 1L << 20;

                auto const freq = 1 / spacing;
                auto const x = uint32_t(offset + freq * pt.x);
                auto const y = uint32_t(offset + freq * pt.y);
                auto const z = uint32_t(offset + freq * pt.z);

                return hash(x, y, z);
            }

        private:
            void init_hash(open_box box)
            {
                // This simple heuristic gives surprisingly good performance.
                hash.modulus = md::linear_hash::uint(box.particle_count * 2 / 11);
                hash.modulus |= 1;
            }

            void init_buckets()
            {
                buckets.resize(hash.modulus);

                // Compute neighbor graph of the spatial cells.
                std::vector<uint32_t> hash_deltas;

                uint32_t const coord_deltas[] = {
                    hash.modulus - 1,
                    hash.modulus,
                    hash.modulus + 1
                };

                for (auto const dx : coord_deltas) {
                    for (auto const dy : coord_deltas) {
                        for (auto const dz : coord_deltas) {
                            hash_deltas.push_back(hash(dx, dy, dz));
                        }
                    }
                }

                sort_unique(hash_deltas);

                for (uint32_t center = 0; center < hash.modulus; center++) {
                    auto& directed_neighbors = buckets[center].directed_neighbors;
                    auto& complete_neighbors = buckets[center].complete_neighbors;

                    for (auto const delta : hash_deltas) {
                        auto const neighbor = (center + delta) % hash.modulus;

                        // Leverage symmetry to reduce search space.
                        if (neighbor >= center) {
                            directed_neighbors.push_back(neighbor);
                        }
                        complete_neighbors.push_back(neighbor);
                    }

                    sort_unique(directed_neighbors);
                    sort_unique(complete_neighbors);
                }
            }
        };
    }
}

#endif
