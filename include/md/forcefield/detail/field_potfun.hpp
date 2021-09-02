// Copyright snsinfu 2018-2021.
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
        template<
            typename P,
            typename = decltype(std::declval<P>()(
                md::index{}
            ))
        >
        void detect_factory_i();

        template<
            typename P,
            typename = decltype(std::declval<P>()(
                std::declval<md::system const&>(),
                md::index{}
            ))
        >
        void detect_factory_si();

        template<typename P, typename = void>
        struct field_potential_factory
        {
            P potential;

            P operator()(md::system const&, md::index) const
            {
                return potential;
            }
        };

        template<typename P>
        struct field_potential_factory<P, decltype(detect_factory_i<P>())>
        {
            P factory;

            auto operator()(md::system const&, md::index i) const
            {
                return factory(i);
            }
        };

        template<typename P>
        struct field_potential_factory<P, decltype(detect_factory_si<P>())>
        {
            P factory;

            auto operator()(md::system const& system, md::index i) const
            {
                return factory(system, i);
            }
        };

        template<typename P>
        field_potential_factory<P> make_field_potential_factory(P pot)
        {
            return field_potential_factory<P>{pot};
        }
    }
}

#endif
