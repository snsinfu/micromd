// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_TRIPLE_POTFUN_HPP
#define MD_FORCEFIELD_DETAIL_TRIPLE_POTFUN_HPP

// This module implements a internal utility to normalize triple potential
// functor into a potential functor factory. Used by make_* family of functions.

#include <utility>

#include "../../basic_types.hpp"


namespace md
{
    namespace detail
    {
        template<typename P>
        struct triple_potential_factory
        {
            P potential;

            inline P operator()(md::index, md::index, md::index) const
            {
                return potential;
            }
        };

        template<
            typename PotFun,
            typename = decltype(std::declval<PotFun>()(md::index{}, md::index{}, md::index{}))
        >
        PotFun make_triple_potfun(PotFun potfun)
        {
            return potfun;
        }

        template<typename P, typename... Dummy>
        triple_potential_factory<P> make_triple_potfun(P pot, Dummy...)
        {
            return triple_potential_factory<P>{pot};
        }
    }
}

#endif
