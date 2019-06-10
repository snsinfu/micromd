#include <iostream>
#include <random>
#include <vector>

#include <md/all.hpp>


int main()
{
    md::system system;
    std::mt19937_64 random;

    std::vector<md::index> reds;
    std::vector<md::index> blacks;

    for (md::index i = 0; i < 5000; i++) {
        std::uniform_real_distribution<md::scalar> normal{0, 0.1};
        system.add_particle(md::basic_particle_data{
            .position = {normal(random), normal(random), normal(random)}
        });

        if (i % 2 == 0) {
            reds.push_back(i);
        } else {
            blacks.push_back(i);
        }
    }

    // Red and black particles repel each other.
    system.add_forcefield(
        md::make_inter_subsystem_neighbor_pair_forcefield(
            md::polybell_potential<8, 8> {
                .overlap_energy  = 10,
                .cutoff_distance = 0.1
            }
        )
        .set_neighbor_distance(0.1)
        .set_subsystems(reds, blacks)
    );

    // Spherical container.
    system.add_forcefield(
        md::make_sphere_outward_forcefield(
            md::harmonic_potential {
                .spring_constant = 1000
            }
        )
        .set_sphere({
            .center = {0, 0, 0},
            .radius = 0.5
        })
    );

    md::simulate_brownian_dynamics(system, {
        .timestep  = 1.0e-6,
        .steps     = 100000,
        .seed      = random(),
        .callback  = [&](md::step step) {
            if (step % 1000 == 0) {
                auto const mean_energy = system.compute_energy() / system.particle_count();
                std::clog << step << '\t' << mean_energy << '\n';
            }
        }
    });

    for (auto part : system.particles()) {
        std::cout << part.position << '\n';
    }
}
