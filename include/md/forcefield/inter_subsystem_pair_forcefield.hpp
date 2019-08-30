// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_FORCEFIELD_INTER_SUBSYSTEM_PAIR_FORCEFIELD_HPP
#define MD_FORCEFIELD_INTER_SUBSYSTEM_PAIR_FORCEFIELD_HPP

// This module provides a template forcefield implementation that computes
// interactions between pairs of particles from two subsystems.

#include <vector>

#include "../basic_types.hpp"
#include "../forcefield.hpp"
#include "../system.hpp"

#include "detail/pair_potfun.hpp"


namespace md
{
    // inter_subsystem_pair_forcefield implements md::forcefield. It computes
    // interactions between pairs of particles in a subsystem.
    //
    // This is a CRTP base class. Derived class must define a callback:
    //
    //     auto inter_subsystem_pair_potential(
    //         md::system const& system,
    //         md::index i,
    //         md::index j
    //     )
    //     Returns the potential object for (i,j) pair.
    //
    template<typename Derived>
    class inter_subsystem_pair_forcefield : public virtual md::forcefield
    {
    public:
        // set_subsystems sets the interacting subsystems.
        Derived& set_subsystems(
            md::array_view<md::index const> sub1,
            md::array_view<md::index const> sub2
        )
        {
            subsystem1_.clear();
            subsystem2_.clear();
            subsystem1_.assign(sub1.begin(), sub1.end());
            subsystem2_.assign(sub2.begin(), sub2.end());
            return derived();
        }

        // compute_energy implements md::forcefield.
        md::scalar compute_energy(md::system const& system) override
        {
            std::vector<particle_eval> part_evals1 = make_part_evals(system, subsystem1_);
            std::vector<particle_eval> part_evals2 = make_part_evals(system, subsystem2_);

            md::scalar sum = 0;

            for (auto& part_i : part_evals1) {
                for (auto& part_j : part_evals2) {
                    auto const pot = derived().inter_subsystem_pair_potential(
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
            std::vector<particle_eval> part_evals1 = make_part_evals(system, subsystem1_);
            std::vector<particle_eval> part_evals2 = make_part_evals(system, subsystem2_);

            for (auto& part_i : part_evals1) {
                for (auto& part_j : part_evals2) {
                    auto const pot = derived().inter_subsystem_pair_potential(
                        system, part_i.index, part_j.index
                    );
                    auto const r = part_i.position - part_j.position;
                    auto const f = pot.evaluate_force(r);
                    part_i.force += f;
                    part_j.force -= f;
                }
            }

            for (particle_eval const& part : part_evals1) {
                forces[part.index] += part.force;
            }
            for (particle_eval const& part : part_evals2) {
                forces[part.index] += part.force;
            }
        }

    private:
        // TODO: Investigate if copying subsystem (like done here) is actually
        // beneficial or not.
        struct particle_eval
        {
            md::index  index    = 0;
            md::point  position = {};
            md::vector force    = {};
        };

        std::vector<particle_eval> make_part_evals(
            md::system const& system, std::vector<md::index> const& targets
        ) const
        {
            auto const positions = system.view_positions();

            std::vector<particle_eval> part_evals;
            part_evals.reserve(targets.size());

            for (md::index idx : targets) {
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

        std::vector<md::index> subsystem1_;
        std::vector<md::index> subsystem2_;
    };

    template<typename PotFun>
    class basic_inter_subsystem_pair_forcefield
        : public md::inter_subsystem_pair_forcefield<
            basic_inter_subsystem_pair_forcefield<PotFun>
        >
    {
    public:
        explicit basic_inter_subsystem_pair_forcefield(PotFun potfun)
            : potfun_{potfun}
        {
        }

        auto inter_subsystem_pair_potential(md::system const&, md::index i, md::index j) const
        {
            return potfun_(i, j);
        }

    private:
        PotFun potfun_;
    };

    // make_inter_subsystem_pair_forcefield implements md::subsystem_pair_forcefield
    // with given potential object or lambda returning a potential object.
    template<typename P>
    auto make_inter_subsystem_pair_forcefield(P pot)
    {
        auto potfun = detail::make_pair_potfun(pot);
        using potfun_type = decltype(potfun);
        return md::basic_inter_subsystem_pair_forcefield<potfun_type>{potfun};
    }
}

#endif
