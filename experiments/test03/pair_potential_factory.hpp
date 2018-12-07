#pragma once

#include <utility>

#include <md/basic_types.hpp>


namespace md
{
    namespace detail
    {
        template<typename P>
        struct pair_potential_factory
        {
            P potential;

            inline P operator()(md::index, md::index) const
            {
                return potential;
            }
        };

        template<
            typename PF,
            typename = decltype(std::declval<PF>()(md::index{}, md::index{}))
        >
        PF make_pair_potfun(PF potfun)
        {
            return potfun;
        }

        template<typename P, typename... Dummy>
        pair_potential_factory<P> make_pair_potfun(P pot, Dummy...)
        {
            return pair_potential_factory<P>{pot};
        }
    }
}
