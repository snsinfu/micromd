// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_FIELD_POTFUN_HPP
#define MD_FORCEFIELD_DETAIL_FIELD_POTFUN_HPP

// This module implements an internal utility to normalize field potential
// functor into a potential functor factory. Used by make_* family of functions.

#include <utility>

#include "../../basic_types.hpp"


namespace md
{
    namespace detail
    {
        template<typename P>
        struct field_potential_factory
        {
            P potential;

            inline P operator()(md::index) const
            {
                return potential;
            }
        };

        template<
            typename PotFun,
            typename = decltype(std::declval<PotFun>()(md::index{}))
        >
        PotFun make_field_potfun(PotFun potfun)
        {
            return potfun;
        }

        template<typename P, typename... Dummy>
        field_potential_factory<P> make_field_potfun(P pot, Dummy...)
        {
            return field_potential_factory<P>{pot};
        }
    }
}

#endif
