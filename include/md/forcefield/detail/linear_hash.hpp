// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_LINEAR_HASH_HPP
#define MD_FORCEFIELD_DETAIL_LINEAR_HASH_HPP

// This module provides a functor class to linearly hash integral triples. Used
// to implement neighbor_searcher.

#include <cstdint>


namespace md
{
    // linear_hash is a linear hash function for integral 3-vectors.
    struct linear_hash
    {
        using hash_t = std::uint32_t;

        hash_t x_coeff = 3929498747;
        hash_t y_coeff = 1008281837;
        hash_t z_coeff = 1832832077;
        hash_t modulus = 1021;

        inline hash_t operator()(hash_t x, hash_t y, hash_t z) const
        {
            // Avoid 32-bit wraparound (mod 2^32) with 64-bit extension.
            using hash2x_t = std::uint64_t;

            hash2x_t sum = 0;
            sum += hash2x_t(x_coeff) * x;
            sum += hash2x_t(y_coeff) * y;
            sum += hash2x_t(z_coeff) * z;

            return hash_t(sum % modulus);
        }
    };
}

#endif
