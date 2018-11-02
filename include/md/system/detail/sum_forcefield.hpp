// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_DETAIL_SUM_FORCEFIELD_HPP
#define MD_SYSTEM_DETAIL_SUM_FORCEFIELD_HPP

#include <cassert>
#include <memory>
#include <vector>

#include "../../basic_types.hpp"
#include "../../forcefield.hpp"


namespace md
{
    class system;

    namespace detail
    {
        // sum_forcefield is a forcefield implementation that compute the sum of
        // other forcefield instances.
        class sum_forcefield : public md::forcefield
        {
        public:
            // add adds a forcefield instance to the sum.
            void add(std::shared_ptr<md::forcefield> ff)
            {
                assert(ff);
                components_.push_back(ff);
            }

            md::scalar compute_energy(md::system const& system) override
            {
                md::scalar sum = 0;

                for (std::shared_ptr<md::forcefield>& component : components_) {
                    sum += component->compute_energy(system);
                }
                return sum;
            }

            void compute_force(md::system const& system, md::array_view<md::vector> forces) override
            {
                for (std::shared_ptr<md::forcefield>& component : components_) {
                    component->compute_force(system, forces);
                }
            }

        private:
            std::vector<std::shared_ptr<md::forcefield>> components_;
        };
    }
}

#endif
