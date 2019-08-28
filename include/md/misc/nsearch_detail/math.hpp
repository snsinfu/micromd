// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_MISC_NSEARCH_DETAIL_MATH_HPP
#define MD_MISC_NSEARCH_DETAIL_MATH_HPP

// This module provides mathematical functions used in neighbor search
// implementations.

#include <cmath>
#include <cstdint>

#include "../../basic_types.hpp"


namespace md
{
    namespace nsearch_detail
    {
        using std::int32_t;
        using std::uint32_t;

        // Computes non-negative floating-point remainder of x/m. The result r
        // satisfies 0 <= r < m. Use this function to map a point to its
        // canonical image in a periodic system.
        inline md::scalar floor_mod(md::scalar x, md::scalar m)
        {
            return x - std::floor(x * (1 / m)) * m;
        }

        // Computes zero-centered floating-point remainder of x/m. The result r
        // satisfies -m/2 <= r <= m/2. Use this function to compute the distance
        // between nearest images of points in a periodic system.
        inline md::scalar round_mod(md::scalar x, md::scalar m)
        {
            return x - std::nearbyint(x * (1 / m)) * m;
        }

        // Returns the integral part of x as uint32_t.
        inline uint32_t trunc_uint(md::scalar x)
        {
            // Converting floating-point number to int32_t is much faster than
            // to uint32_t on x86.
            return static_cast<uint32_t>(static_cast<int32_t>(x));
        }

        // Rounds x to the nearest uint32_t value.
        inline uint32_t round_uint(md::scalar x)
        {
            return trunc_uint(std::nearbyint(x));
        }

        // Quickly computes (x + y) mod m assuming both x and y are less than m.
        inline uint32_t add_mod(uint32_t x, uint32_t y, uint32_t m)
        {
            auto sum = x + y;
            if (sum >= m) {
                sum -= m;
            }
            return sum;
        }

        // A bin_layout defines 1-dimensional uniform bins.
        struct bin_layout
        {
            md::scalar step;
            uint32_t count;
        };

        // Creates a bin_layout that splits given span into uniform bins.
        inline bin_layout define_bins(md::scalar span, md::scalar min_step)
        {
            bin_layout bins;
            bins.count = std::max(trunc_uint(span / min_step), uint32_t(1));
            bins.step = span / md::scalar(bins.count);
            return bins;
        }

        // Computes the index of the bin in which given coordinate value falls.
        inline uint32_t locate_bin(const bin_layout& bins, md::scalar coord)
        {
            auto const span = bins.step * md::scalar(bins.count);
            auto const image = floor_mod(coord, span);
            auto const index = trunc_uint(image * (1 / bins.step));
            return index < bins.count ? index : bins.count - 1;
        }
    }
}

#endif
