// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_SUBSYSTEM_PAIR_FORCEFIELD_HPP
#define MD_FORCEFIELD_SUBSYSTEM_PAIR_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// interactions between particle pairs in a subsystem.

#include <type_traits>
#include <utility>
#include <vector>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/pair_potfun.hpp"


namespace md
{
    // subsystem_pair_forcefield implements md::forcefield. It computes
    // interactions between pairs of particles in a subsystem.
    //
    // This is a CRTP base class. Derived class must define a callback:
    //
    //     auto subsystem_pair_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) pair.
    //
    template<typename Derived>
    class subsystem_pair_forcefield : public virtual md::forcefield
    {
    public:
        // add_target selects a particle.
        Derived& add_target(md::index idx)
        {
            targets_.push_back(idx);
            return derived();
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            std::vector<particle_eval> part_evals = make_part_evals(system);
            md::scalar sum = 0;

            for (md::index i = 0; i < part_evals.size(); i++) {
                for (md::index j = i + 1; j < part_evals.size(); j++) {
                    auto& part_i = part_evals[i];
                    auto& part_j = part_evals[j];

                    auto const pot = derived().subsystem_pair_potential(
                        system, part_i.index, part_j.index
                    );
                    auto const r = part_i.position - part_j.position;

                    sum += pot.evaluate_energy(r);
                }
            }

            return sum;
        }

        // compute_force implements md::forcefield.
        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            std::vector<particle_eval> part_evals = make_part_evals(system);

            for (md::index i = 0; i < part_evals.size(); i++) {
                for (md::index j = i + 1; j < part_evals.size(); j++) {
                    auto& part_i = part_evals[i];
                    auto& part_j = part_evals[j];

                    auto const pot = derived().subsystem_pair_potential(
                        system, part_i.index, part_j.index
                    );
                    auto const r = part_i.position - part_j.position;
                    auto const f = pot.evaluate_force(r);
                    part_i.force += f;
                    part_j.force -= f;
                }
            }

            for (particle_eval const& part : part_evals) {
                forces[part.index] += part.force;
            }
        }

    private:
        // Double indirections like positions[targets_[i]] harm performance. So
        // copy subsystem into a temp array and work on it.
        struct alignas(64) particle_eval
        {
            md::index  index    = 0;
            md::point  position = {};
            md::vector force    = {};
        };

        std::vector<particle_eval> make_part_evals(md::system const& system) const
        {
            auto const positions = system.view_positions();

            std::vector<particle_eval> part_evals;
            part_evals.reserve(targets_.size());

            for (md::index idx : targets_) {
                particle_eval part;
                part.index = idx;
                part.position = positions[idx];
                part_evals.push_back(part);
            }

            return part_evals;
        }

        // derived returns a reference to this object as the CRTP derived class.
        Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }

        std::vector<md::index> targets_;
    };

    template<typename PotFun>
    class basic_subsystem_pair_forcefield
        : public md::subsystem_pair_forcefield<basic_subsystem_pair_forcefield<PotFun>>
    {
    public:
        explicit basic_subsystem_pair_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto subsystem_pair_potential(md::system const&, md::index i, md::index j) const
        {
            return potfun_(i, j);
        }

    private:
        PotFun potfun_;
    };

    // make_subsystem_pair_forcefield implements md::subsystem_pair_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename P>
    auto make_subsystem_pair_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_subsystem_pair_forcefield<potfun_type>{potfun};
    }
}

#endif
