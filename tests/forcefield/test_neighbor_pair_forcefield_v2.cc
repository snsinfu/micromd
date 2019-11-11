#include <algorithm>
#include <cassert>
#include <iterator>
#include <set>
#include <utility>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>
#include <md/potential/softcore_potential.hpp>

#include <md/forcefield/neighbor_pair_forcefield_v2.hpp>

#include <catch.hpp>


namespace
{
    md::scalar max_difference(
        md::array_view<md::vector const> v1,
        md::array_view<md::vector const> v2
    )
    {
        assert(v1.size() == v2.size());

        md::scalar diff = 0;

        for (md::index i = 0; i < v1.size(); i++) {
            diff = std::max(diff, md::norm(v1[i] - v2[i]));
        }
        return diff;
    }
}

TEST_CASE("neighbor_pair_forcefield_v2 - computes correct forcefield")
{
    md::scalar const cutoff_distance = 0.3;
    md::index const point_count = 1000;

    md::periodic_box box;
    box.x_period = 0.9;
    box.y_period = 1.0;
    box.z_period = 1.1;

    // Random points in a box. Coordinate values are deliberately overdispersed
    // compared to the periods to test handling of periodic boundary conditions.
    md::system system;
    std::mt19937 random;
    std::uniform_real_distribution<md::scalar> coord{-3, 3};
    for (md::index i = 0; i < point_count; i++) {
        auto part = system.add_particle();
        part.position = { coord(random), coord(random), coord(random) };
    }

    // Short-range interactions.
    md::softcore_potential<3> potential;
    potential.overlap_energy = 1.0;
    potential.cutoff_distance = cutoff_distance;

    auto forcefield = md::make_neighbor_pair_forcefield_v2<md::periodic_box>(
        [&](md::index, md::index) {
            return potential;
        }
    )
    .set_box(box)
    .set_neighbor_distance(cutoff_distance);

    // Test forcefield.
    std::vector<md::vector> actual_forces(system.particle_count());
    std::vector<md::vector> expect_forces(system.particle_count());
    forcefield.compute_force(system, actual_forces);

    md::scalar actual_energy = forcefield.compute_energy(system);
    md::scalar expect_energy = 0;

    // Ground truth by brute-force loop.
    md::array_view<md::point const> positions = system.view_positions();

    for (md::index i = 0; i < positions.size(); i++) {
        for (md::index j = i + 1; j < positions.size(); j++) {
            md::vector const r = box.shortest_displacement(positions[i], positions[j]);
            md::vector const force = potential.evaluate_force(r);
            md::scalar const energy = potential.evaluate_energy(r);
            expect_forces[i] += force;
            expect_forces[j] -= force;
            expect_energy += energy;
        }
    }

    CHECK(actual_energy == Approx(expect_energy));
    CHECK(max_difference(actual_forces, expect_forces) == Approx(0).margin(0.001));
}

TEST_CASE("neighbor_pair_forcefield_v2::set_targets - limits search targets")
{
    md::scalar const cutoff_distance = 0.1;

    md::periodic_box box;
    box.x_period = 0.9;
    box.y_period = 1.0;
    box.z_period = 1.1;

    std::set<md::index> const targets = {
        2, 3, 5, 7
    };

    // Most particles are within the cutoff distance.
    md::system system;

    system.add_particle().position = {0.1, 0.1, 0.1};
    system.add_particle().position = {0.1, 0.1, 0.1};
    system.add_particle().position = {0.1, 0.1, 0.1}; // 2
    system.add_particle().position = {0.1, 0.1, 0.1}; // 3
    system.add_particle().position = {0.1, 0.1, 0.1};
    system.add_particle().position = {0.0, 0.0, 0.0}; // 5
    system.add_particle().position = {0.0, 0.0, 0.0};
    system.add_particle().position = {0.0, 0.0, 0.0}; // 7
    system.add_particle().position = {0.0, 0.0, 0.0};
    system.add_particle().position = {0.0, 0.0, 0.0};

    // Short-range interactions.
    md::softcore_potential<3> potential;
    potential.overlap_energy = 1.0;
    potential.cutoff_distance = cutoff_distance;

    auto forcefield = md::make_neighbor_pair_forcefield_v2<md::periodic_box>(
        potential
    )
    .set_targets(targets)
    .set_box(box)
    .set_neighbor_distance(cutoff_distance);

    md::array_view<md::point const> positions = system.view_positions();
    md::scalar actual_energy = forcefield.compute_energy(system);
    md::scalar expect_energy = 0;

    for (md::index const i : targets) {
        for (md::index const j : targets) {
            if (i >= j) {
                continue;
            }
            md::vector const r = box.shortest_displacement(positions[i], positions[j]);
            if (r.norm() < cutoff_distance) {
                expect_energy += potential.evaluate_energy(r);
            }
        }
    }

    CHECK(actual_energy == Approx(expect_energy));
}
