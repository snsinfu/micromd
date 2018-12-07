#pragma once

#include <md/basic_types.hpp>
#include <md/system.hpp>

#include <md/forcefield/sequential_pair_forcefield.hpp>


namespace md
{
    template<typename PotFun>
    class basic_sequential_pair_forcefield
        : public md::sequential_pair_forcefield<basic_sequential_pair_forcefield<PotFun>>
    {
    public:
        explicit basic_sequential_pair_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto sequential_pair_potential(md::system const&, md::index i, md::index j) const
        {
            return potfun_(i, j);
        }

    private:
        md::scalar neighbor_distance_ = 0;
        PotFun potfun_;
    };

    template<typename P>
    auto make_sequential_pair_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        return basic_sequential_pair_forcefield<decltype(potfun)>{potfun};
    }
}
