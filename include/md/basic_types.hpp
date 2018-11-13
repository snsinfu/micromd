// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_BASIC_TYPES_HPP
#define MD_BASIC_TYPES_HPP

// This module defines basic typedefs used across micromd library.

#include <cstddef>
#include <cstdint>

#include "basic_types/array_view.hpp"
#include "basic_types/point.hpp"


namespace md
{
    // index is the integral type used to index arrays.
    using index = std::size_t;

    // scalar is the floating-point type of choice.
    using scalar = double;

    // step is the integral type used to count simulation steps.
    using step = std::int64_t;
}

#endif
