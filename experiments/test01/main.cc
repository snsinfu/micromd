#include <iostream>
#include <vector>

#include <md/all.hpp>


class particle_forcefield : public md::point_source_forcefield<particle_forcefield>
{
public:
    auto point_source_potential(md::system const&, md::index) const
    {
        return md::harmonic_potential{1.234};
    }
};


int main()
{
    md::system system;

    system.add_particle({
        .mobility = 1.0,
        .position = {0, 0, 0},
    });

    system.add_forcefield(particle_forcefield{});

    auto callback = [&](md::step step) {
        std::cout << system.view_positions()[0] << '\n';
    };

    md::simulate_brownian_dynamics(system, {
        .temperature = 1,
        .timestep    = 1,
        .steps       = 5000,
        .callback    = callback,
    });
}
