#pragma once

#include <md/basic_types.hpp>
#include <md/system.hpp>

#include <md/forcefield/neighbor_pair_forcefield.hpp>


namespace md
{
    template<typename PotFun>
    class basic_neighbor_pair_forcefield
        : public md::neighbor_pair_forcefield<basic_neighbor_pair_forcefield<PotFun>>
    {
    public:
        explicit basic_neighbor_pair_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        basic_neighbor_pair_forcefield& set_neighbor_distance(md::scalar dist)
        {
            neighbor_distance_ = dist;
            return *this;
        }

        md::scalar neighbor_distance(md::system const&) const
        {
            return neighbor_distance_;
        }

        auto neighbor_pair_potential(md::system const&, md::index i, md::index j) const
        {
            return potfun_(i, j);
        }

    private:
        md::scalar neighbor_distance_ = 0;
        PotFun potfun_;
    };

    template<typename P>
    auto make_neighbor_pair_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        return basic_neighbor_pair_forcefield<decltype(potfun)>{potfun};
    }
}
