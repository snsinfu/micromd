#include <iostream>
#include <vector>

#include <md/all.hpp>


class particle_forcefield : public md::composite_forcefield<
    md::neighbor_pair_forcefield  <particle_forcefield>,
    md::sequential_pair_forcefield<particle_forcefield>
>
{
public:
    auto neighbor_distance(md::system const&) const
    {
        return 0.3;
    }

    auto neighbor_pair_potential(md::system const&, md::index, md::index) const
    {
        return md::softcore_potential<4>{
            .overlap_energy  = 4.0,
            .cutoff_distance = 0.3,
        };
    }

    auto sequential_pair_potential(md::system const&, md::index, md::index) const
    {
        return md::harmonic_potential{
            .spring_constant = 25,
        };
    }
};


int main()
{
    md::system system;
    particle_forcefield forcefield;

    for (int i = 0; i < 10; i++) {
        system.add_particle({
            .position = {0.3 * i, 0, 0}
        });
    }

    forcefield.add_segment(0, 9);

    system.add_forcefield(forcefield);

    auto callback = [&](md::step step) {
        if ((step + 1) % 1000 == 0) {
            std::cout
                << system.view_positions()[1] - system.view_positions()[0]
                << '\n';
        }
    };

    md::simulate_brownian_dynamics(system, {
        .temperature = 1,
        .timestep    = 0.0026,
        .steps       = 500000,
        .callback    = callback,
    });
}
