#include <iostream>
#include <vector>

#include <md/all.hpp>

#include "make_neighbor_pair_forcefield.hpp"
#include "make_sequential_pair_forcefield.hpp"


int main()
{
    md::system system;

    for (int i = 0; i < 10; i++) {
        system.add_particle({
            .position = {0.3 * i, 0, 0}
        });
    }

    system.add_forcefield(
        md::make_neighbor_pair_forcefield(md::softcore_potential<4>{
            .overlap_energy  = 4.0,
            .cutoff_distance = 0.3,
        })
        .set_neighbor_distance(0.3)
    );

    system.add_forcefield(
        md::make_sequential_pair_forcefield(md::harmonic_potential{
            .spring_constant = 40,
        })
    );

    auto callback = [&](md::step step) {
        if ((step + 1) % 1000 == 0) {
            std::cout
                << system.view_positions()[1] - system.view_positions()[0]
                << '\n';
        }
    };

    md::simulate_brownian_dynamics(system, {
        .temperature = 1,
        .timestep    = 1e-5,
        .steps       = 500000,
        .callback    = callback,
    });
}
