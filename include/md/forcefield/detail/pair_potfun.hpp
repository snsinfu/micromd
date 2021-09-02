// Copyright snsinfu 2018-2021.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_DETAIL_PAIR_POTFUN_HPP
#define MD_FORCEFIELD_DETAIL_PAIR_POTFUN_HPP

// This module implements a internal utility to normalize pair potential or
// pair potential functor into a potential factory object. Used by make_*
// family of forcefield functions.

#include <utility>

#include "../../basic_types.hpp"


namespace md
{
    namespace detail
    {
        template<
            typename P,
            typename = decltype(std::declval<P>()(
                md::index{},
                md::index{}
            ))
        >
        void detect_factory_ij();

        template<
            typename P,
            typename = decltype(std::declval<P>()(
                std::declval<md::system const&>(),
                md::index{},
                md::index{}
            ))
        >
        void detect_factory_sij();

        template<typename P, typename = void>
        struct pair_potential_factory
        {
            P potential;

            P operator()(md::system const&, md::index, md::index) const
            {
                return potential;
            }
        };

        template<typename P>
        struct pair_potential_factory<P, decltype(detect_factory_ij<P>())>
        {
            P factory;

            auto operator()(md::system const&, md::index i, md::index j) const
            {
                return factory(i, j);
            }
        };

        template<typename P>
        struct pair_potential_factory<P, decltype(detect_factory_sij<P>())>
        {
            P factory;

            auto operator()(md::system const& system, md::index i, md::index j) const
            {
                return factory(system, i, j);
            }
        };

        template<typename P>
        pair_potential_factory<P> make_pair_potential_factory(P pot)
        {
            return pair_potential_factory<P>{pot};
        }
    }
}

#endif
