#include <cmath>
#include <numeric>
#include <vector>

#include <md/basic_types.hpp>
#include <md/system.hpp>
#include <md/forcefield/all_pair_forcefield.hpp>
#include <md/forcefield/sequential_pair_forcefield.hpp>
#include <md/forcefield/sequential_triple_forcefield.hpp>
#include <md/potential/cosine_bending_potential.hpp>
#include <md/potential/softcore_potential.hpp>
#include <md/potential/spring_potential.hpp>
#include <md/simulation/brownian_dynamics.hpp>

#include <catch.hpp>

#include <iostream>


namespace
{
    md::scalar compute_slope(
        std::vector<md::scalar> const& xs,
        std::vector<md::scalar> const& ys
    )
    {
        md::scalar mean_x = std::accumulate(xs.begin(), xs.end(), md::scalar(0));
        mean_x /= md::scalar(xs.size());

        md::scalar mean_y = std::accumulate(ys.begin(), ys.end(), md::scalar(0));
        mean_y /= md::scalar(ys.size());

        md::scalar xx = 0;
        md::scalar xy = 0;

        for (md::index i = 0; i < xs.size(); i++) {
            xx += (xs[i] - mean_x) * (xs[i] - mean_x);
            xy += (xs[i] - mean_x) * (ys[i] - mean_y);
        }

        return xy / xx;
    }
}


TEST_CASE("Persistence length of a polymer", "[.][slow]")
{
    // Consider a homopolymer. Let bi be the i-th bond vector. The canonical
    // mean of the cosine of the angle between two bond vectors bi and bj
    // decays exponentially with respect to the contour distance:
    //
    //     E[cos(bi, bj)] = exp(-L/P) ,   L = |i - j| .
    //
    // The characteristic decay distance P is called the persistence length of
    // the polymer. If bending energy of the polymer is B, the persistence
    // length would be P = B / kT.

    md::scalar const B = 20;
    md::scalar const kT = 1;
    md::scalar const a = 0.1;
    md::index const max_L = 10;

    md::system system;

    for (md::index i = 0; i < 2 * max_L; i++) {
        system.add_particle().position = {a * i, 0, 0};
    }

    system.add_forcefield(
        md::make_all_pair_forcefield(
            md::softcore_potential<4>{50, a}
        )
    );

    system.add_forcefield(
        md::make_sequential_pair_forcefield(
            md::spring_potential{400, a}
        )
        .add_segment(0, system.particle_count() - 1)
    );

    system.add_forcefield(
        md::make_sequential_triple_forcefield(
            md::cosine_bending_potential{B}
        )
        .add_segment(0, system.particle_count() - 1)
    );

    // Sampling
    md::brownian_dynamics_config config;
    config.temperature = kT;
    config.timestep = 1e-4;
    config.steps = 100000;

    std::vector<md::scalar> mean_cosines(max_L);

    config.callback = [&](md::step step) {
        md::array_view<md::point const> const pos = system.view_positions();
        md::vector const b1 = pos[1] - pos[0];

        for (md::index L = 0; L < max_L; L++) {
            md::vector const bi = pos[L + 1] - pos[L];
            md::scalar const cos = md::dot(b1.normalize(), bi.normalize());

            mean_cosines[L] += (cos - mean_cosines[L]) / md::scalar(step);
        }
    };

    md::simulate_brownian_dynamics(system, config);

    // Fit y = a*x to the data: (L, log(mean_cos)).
    std::vector<md::scalar> xs;
    std::vector<md::scalar> ys;

    for (md::index L = 0; L < max_L; L++) {
        xs.push_back(md::scalar(L));
        ys.push_back(std::log(mean_cosines[L]));
    }

    md::scalar const slope = compute_slope(xs, ys);
    md::scalar const P = -1 / slope;

    CHECK(P == Approx(B / kT).epsilon(0.1));
}
