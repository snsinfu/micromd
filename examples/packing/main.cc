#include <chrono>
#include <iostream>
#include <random>

#include <md/all.hpp>


namespace
{
    md::scalar implicit(md::sphere const& sphere, md::point const& pt)
    {
        return (pt - sphere.center).squared_norm() - sphere.radius * sphere.radius;
    }

    md::scalar implicit(md::ellipsoid const& el, md::point const& pt)
    {
        md::vector const r = pt - el.center;
        md::vector const e = {
            1 / (el.semiaxis_x * el.semiaxis_x),
            1 / (el.semiaxis_y * el.semiaxis_y),
            1 / (el.semiaxis_z * el.semiaxis_z)
        };
        return dot(e.hadamard(r), r) - 1;
    }

    class packing_forcefield : public md::composite_forcefield<
        md::neighbor_pair_forcefield<packing_forcefield>,
        md::ellipsoid_surface_forcefield<packing_forcefield>
    >
    {
    public:
        auto neighbor_distance(md::system const&) const
        {
            return 0.1;
        }

        auto neighbor_pair_potential(md::system const&, md::index, md::index) const
        {
            md::softcore_potential<4> pot;
            pot.overlap_energy = 6.0;
            pot.cutoff_distance = 0.1;
            return pot;
        }

        auto ellipsoid_outward_potential(md::system const&, md::index) const
        {
            md::harmonic_potential pot;
            pot.spring_constant = 2000;
            return pot;
        }
    };
}


int main()
{
    md::system system;

    std::mt19937_64 random;
    std::uniform_real_distribution<md::scalar> coord{-4, 4};

    md::ellipsoid ellipsoid;
    ellipsoid.semiaxis_x = 4;
    ellipsoid.semiaxis_y = 1;
    ellipsoid.semiaxis_z = 0.25;

    for (int i = 0; i < 1000; i++) {
        auto particle = system.add_particle();
        do {
            particle.position = {
                coord(random),
                coord(random),
                coord(random)
            };
        } while (implicit(ellipsoid, particle.position) > 0);
    }

    packing_forcefield forcefield;
    forcefield.set_ellipsoid(ellipsoid);
    system.add_forcefield(forcefield);

    md::brownian_dynamics_config config;

    config.temperature = 1;
    config.timestep = 1e-4;
    config.spacestep = 0.01;
    config.steps = 500000;

    using clock = std::chrono::steady_clock;
    auto logging_interval = md::step(1000);
    auto prev_timestamp = clock::now();

    config.callback = [&](md::step step) {
        if ((step + 1) % logging_interval == 0) {
            auto timestamp = clock::now();
            auto time = std::chrono::duration_cast<std::chrono::duration<double>>(
                timestamp - prev_timestamp
            );
            std::clog
                << step + 1
                << '\t'
                << time.count() * 1e3 / logging_interval
                << " ms/step"
                << '\t'
                << system.compute_potential_energy()
                << '\n';
            prev_timestamp = timestamp;
        }
    };

    md::simulate_brownian_dynamics(system, config);
}
