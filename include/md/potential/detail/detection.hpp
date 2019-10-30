// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_DETAIL_DETECTION_HPP
#define MD_POTENTIAL_DETAIL_DETECTION_HPP

#include "../../basic_types.hpp"

namespace md
{
    namespace detail
    {
        // SFINAE utility that detects pairwise potential functor.
        template<
            typename Pot,
            md::scalar(Pot::*)(md::vector) const = &Pot::evaluate_energy,
            md::vector(Pot::*)(md::vector) const = &Pot::evaluate_force
        >
        Pot detect_pairwise_potential();

        // `pairwise_potential_t` aliases the given type `Pot` if `Pot` is a
        // pairwise potential functor. Otherwise instantiation fails.
        template<typename Pot>
        using pairwise_potential_t = decltype(detect_pairwise_potential<Pot>());
    }
}

#endif
