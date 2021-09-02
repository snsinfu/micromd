#include <algorithm>
#include <iostream>
#include <random>

#include <md.hpp>


// Pack 10000 particles in a spherical container.

int main()
{
    md::system system;

    std::mt19937_64 random;
    std::uniform_real_distribution<md::scalar> uniform{-0.5, 0.5};

    for (int i = 0; i < 10000; i++) {
        system.add_particle(md::basic_particle_data{
            .position = {uniform(random), uniform(random), uniform(random)}
        });
    }

    system.add_forcefield(
        md::make_neighbor_pairwise_forcefield(
            md::softcore_potential<2, 3> {
                .energy = 5.0,
                .diameter = 0.1
            }
        )
        .set_neighbor_distance(0.1)
    );

    system.add_forcefield(
        md::make_sphere_outward_forcefield(
            md::harmonic_potential {
                .spring_constant = 1000
            }
        )
        .set_sphere(md::sphere {
            .center = {0, 0, 0},
            .radius = 1
        })
    );

    std::clog << "step\tenergy\n";

    md::simulate_brownian_dynamics(system, {
        .timestep = 1e-5,
        .steps = 10000,
        .callback = [&](md::step step) {
            if (step % 100 == 0) {
                auto const energy =
                system.compute_energy() / md::scalar(system.particle_count());
                std::clog << step << '\t' << energy << '\n';
            }
        }
    });
}
