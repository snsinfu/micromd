// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_TYPEDEF_HPP
#define MD_TYPEDEF_HPP

#include <cstddef>

#include "vendor/array_view.hpp"
#include "vendor/point.hpp"


namespace md
{
    // index is the integral type used to index arrays.
    using index = std::size_t;

    // scalar is the floating-point type of choice.
    using scalar = double;
}

#endif
