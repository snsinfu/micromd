#include <algorithm>
#include <iostream>
#include <random>

#include <md.hpp>


class particle_forcefield : public md::composite_forcefield<
    md::neighbor_pair_forcefield<md::open_box, particle_forcefield>,
    md::sphere_surface_forcefield<particle_forcefield>
>
{
public:
    md::scalar neighbor_distance(md::system const&) const
    {
        return 0.1;
    }

    auto neighbor_pair_potential(md::system const&, md::index, md::index) const
    {
        return md::softcore_potential<2, 3> {
            .energy = 5.0,
            .diameter = 0.1
        };
    }

    auto sphere_outward_potential(md::system const&, md::index) const
    {
        return md::harmonic_potential {
            .spring_constant = 1000,
        };
    }
};


int main()
{
    md::system system;

    std::mt19937_64 random;
    std::uniform_real_distribution<md::scalar> uniform{-1, 1};

    for (int i = 0; i < 10000; i++) {
        system.add_particle(md::basic_particle_data{
            .position = {uniform(random), uniform(random), uniform(random)}
        });
    }

    particle_forcefield forcefield;
    forcefield.set_sphere({
        .center = {0, 0, 0},
        .radius = 1
    });
    system.add_forcefield(forcefield);

    md::simulate_brownian_dynamics(system, {
        .timestep = 1.0e-5,
        .spacestep = 0.001,
        .steps = 10000
    });
}
