#include <iostream>
#include <vector>

#include <md/all.hpp>


struct spring_potential
{
    md::scalar equilibrium_distance = 1;
    md::scalar spring_constant = 0;

    md::scalar evaluate_energy(md::vector r) const
    {
        md::scalar const u = r.norm() - equilibrium_distance;
        return 0.5 * spring_constant * u * u;
    }

    md::vector evaluate_force(md::vector r) const
    {
        md::scalar const d = r.norm();
        if (d == 0) {
            return {};
        }
        return spring_constant * (equilibrium_distance / d - 1) * r;
    }
};


int main()
{
    md::system system;

    for (int i = 0; i < 10; i++) {
        system.add_particle({
            .position = {0.3 * i, 0, 0}
        });
    }

    /*
    system.add_forcefield(
        md::make_neighbor_pair_forcefield(md::softcore_potential<4>{
            .overlap_energy  = 4.0,
            .cutoff_distance = 0.3,
        })
        .set_neighbor_distance(0.3)
    );
    */

    system.add_forcefield(
        md::make_sequential_pair_forcefield(spring_potential{
            .equilibrium_distance = 0.3,
            .spring_constant      = 200,
        })
        .add_segment(0, 9)
    );

    system.add_forcefield(
        md::make_sequential_triple_forcefield(md::cosine_bending_potential{
            .bending_energy = 100,
        })
        .add_segment(0, 9)
    );

    auto callback = [&](md::step step) {
        if (step % 100 == 0) {
            auto const pos = system.view_positions();
            auto const b1 = pos[1] - pos[0];
            auto const b2 = pos[pos.size() - 1] - pos[pos.size() - 2];
            std::cout
                << b1.norm() << '\t'
                << b2.norm() << '\t'
                << dot(b1.normalize(), b2.normalize()) << '\n';
        }
    };

    md::simulate_brownian_dynamics(system, {
        .temperature = 1,
        .timestep    = 1e-6,
        .steps       = 1000000,
        .callback    = callback,
    });
}
