#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <vector>

#include <md/basic_types.hpp>
#include <md/system.hpp>
#include <md/simulation/brownian_dynamics.hpp>

#include <catch.hpp>


TEST_CASE("Sampling MSD of Brownian particle")
{
    md::scalar const kT = 1.23;
    md::scalar const mu = 4.56;
    md::scalar const dt = 1e-4;

    // Sample MSD
    md::system system;

    for (int i = 0; i < 1000; i++) {
        system.add_particle().mobility = mu;
    }

    md::brownian_dynamics_config config;

    config.temperature = kT;
    config.timestep = dt;
    config.steps = 100000;

    std::vector<md::scalar> r2_samples;

    config.callback = [&](md::step step) {
        md::scalar mean_r2 = 0;

        for (md::point const pt : system.view_positions()) {
            mean_r2 += pt.vector().squared_norm();
        }
        mean_r2 /= md::scalar(system.particle_count());

        r2_samples.push_back(mean_r2);
    };

    md::simulate_brownian_dynamics(system, config);

    md::scalar const D = 6 * mu * kT;
}
