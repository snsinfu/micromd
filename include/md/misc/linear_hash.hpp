// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_MISC_LINEAR_HASH_HPP
#define MD_MISC_LINEAR_HASH_HPP

// This module provides a functor class to linearly hash integral pairs. Used
// to implement neighbor_searcher.

#include <cstdint>


namespace md
{
    // linear_hash is a linear hash function for integral pairs. The following
    // equality holds for any linear_hash h:
    //
    //     h(x+dx, y+dy) = h(x,y) + h(dx,dy) .
    //
    struct linear_hash
    {
        using uint = std::uint32_t;

        // Hash coefficients. Default values are arbitrarily chosen primes.
        uint x_coeff = 3929498747;
        uint y_coeff = 1008281837;
        uint modulus = 1021;

        inline uint operator()(uint x, uint y) const
        {
            // Avoid 32-bit wraparound, which breaks linearity.
            using uint2x = std::uint64_t;

            uint2x sum = 0;
            sum += uint2x{x_coeff} * x;
            sum += uint2x{y_coeff} * y;

            return uint(sum % modulus);
        }
    };
}

#endif
