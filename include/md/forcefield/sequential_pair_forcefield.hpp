// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_SEQUENTIAL_PAIR_FORCEFIELD_HPP
#define MD_FORCEFIELD_SEQUENTIAL_PAIR_FORCEFIELD_HPP

#include <utility>
#include <vector>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"


namespace md
{
    template<typename Derived>
    class sequential_pair_forcefield : public virtual md::forcefield
    {
    public:
        void add_segment(md::index first, md::index last)
        {
            segments_.emplace_back(first, last);
        }

        md::scalar compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            md::scalar sum = 0;

            for (std::pair<md::index, md::index> segment : segments_) {
                for (md::index i = segment.first; i < segment.second; i++) {
                    md::index const j = i + 1;
                    md::vector const r = positions[i] - positions[j];

                    auto const pot = derived().sequential_pair_potential(system, i, j);
                    md::scalar const energy = pot.evaluate_energy(r);

                    sum += energy;
                }
            }

            return sum;
        }

        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (std::pair<md::index, md::index> segment : segments_) {
                for (md::index i = segment.first; i < segment.second; i++) {
                    md::index const j = i + 1;
                    md::vector const r = positions[i] - positions[j];

                    auto const pot = derived().sequential_pair_potential(system, i, j);
                    md::vector const force = pot.evaluate_force(r);

                    forces[i] += force;
                    forces[j] -= force;
                }
            }
        }

    private:
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        std::vector<std::pair<md::index, md::index>> segments_;
    };
}

#endif
