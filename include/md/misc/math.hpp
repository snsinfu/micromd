// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_MISC_MATH_HPP
#define MD_MISC_MATH_HPP

// This module provides a bell-shaped polynomial potential.

#include <cmath>

#include "../basic_types.hpp"


namespace md
{
    // Returns x raised to the power of N.
    template<int N>
    md::scalar power(md::scalar x)
    {
        if (N == 0) {
            return 1;
        }

        // Exponentiation by squaring. The loop would be unrolled by compiler.
        md::scalar pow = x;
        for (int n = N - 1; n > 0; n /= 2) {
            if (n % 2) {
                pow *= x;
            }
            x *= x;
        }
        return pow;
    }

    // Returns x raised to the power of N/2.
    template<int N>
    md::scalar power_sqrt(md::scalar x)
    {
        if (N % 2 == 0) {
            return power<N / 2>(x);
        }
        return power<N>(std::sqrt(x));
    }
}

#endif
