// Copyright snsinfu 2018-2021.
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
        template<
            typename P,
            typename = decltype(std::declval<P>()(
                md::index{},
                md::index{},
                md::index{}
            ))
        >
        void detect_factory_ijk();

        template<
            typename P,
            typename = decltype(std::declval<P>()(
                std::declval<md::system const&>(),
                md::index{},
                md::index{},
                md::index{}
            ))
        >
        void detect_factory_sijk();

        template<typename P, typename = void>
        struct triple_potential_factory
        {
            P potential;

            inline
            P operator()(md::system const&, md::index, md::index, md::index) const
            {
                return potential;
            }
        };

        template<typename P>
        struct triple_potential_factory<P, decltype(detect_factory_ijk<P>())>
        {
            P factory;

            inline
            auto operator()(md::system const&, md::index i, md::index j, md::index k) const
            {
                return factory(i, j, k);
            }
        };

        template<typename P>
        struct triple_potential_factory<P, decltype(detect_factory_sijk<P>())>
        {
            P factory;

            inline
            auto operator()(md::system const& system, md::index i, md::index j, md::index k) const
            {
                return factory(system, i, j, k);
            }
        };

        template<typename P>
        triple_potential_factory<P> make_triple_potential_factory(P pot)
        {
            return triple_potential_factory<P>{pot};
        }
    }
}

#endif
