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

#include <md/forcefield/intra_subsystem_neighbor_pair_forcefield.hpp>

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

TEST_CASE("intra_subsystem_neighbor_pair_forcefield - computes correct forcefield")
{
    class test_forcefield : public md::intra_subsystem_neighbor_pair_forcefield<test_forcefield>
    {
    public:
        md::scalar intra_subsystem_neighbor_distance(md::system const&)
        {
            return 1;
        }

        md::softcore_potential<2, 2> intra_subsystem_neighbor_pair_potential(
            md::system const&, md::index, md::index
        )
        {
            return md::softcore_potential<2, 2>{1, 1};
        }
    };

    // Particles on a 5x5x5 grid
    md::system system;
    std::vector<md::index> reds;

    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            for (int z = -2; z <= 2; z++) {
                md::basic_particle_data data;

                data.position = {
                    x * 0.5,
                    y * 0.5,
                    z * 0.5
                };

                md::particle_ref part = system.add_particle(data);

                if (part.index % 2 == 0) {
                    reds.push_back(part.index);
                }
            }
        }
    }

    test_forcefield forcefield;
    forcefield.set_subsystem(reds);

    // Naive calculation.
    md::softcore_potential<2, 2> potential;

    std::vector<md::vector> actual_forces(system.particle_count());
    std::vector<md::vector> expected_forces(system.particle_count());

    forcefield.compute_force(system, actual_forces);

    md::scalar actual_energy = forcefield.compute_energy(system);
    md::scalar expected_energy = 0;

    md::array_view<md::point const> positions = system.view_positions();

    for (md::index i : reds) {
        for (md::index j : reds) {
            if (j <= i) {
                continue;
            }
            md::vector const r = positions[i] - positions[j];
            md::vector const force = potential.evaluate_force(r);
            md::scalar const energy = potential.evaluate_energy(r);

            expected_forces[i] += force;
            expected_forces[j] -= force;

            expected_energy += energy;
        }
    }

    CHECK(actual_energy == Approx(expected_energy));
    CHECK(max_difference(actual_forces, expected_forces) < 1e-6);
}

TEST_CASE("intra_subsystem_neighbor_pair_forcefield::compute_force - adds force to array")
{
    class test_forcefield : public md::intra_subsystem_neighbor_pair_forcefield<test_forcefield>
    {
    public:
        md::scalar intra_subsystem_neighbor_distance(md::system const&)
        {
            return 1;
        }

        md::softcore_potential<2, 2> intra_subsystem_neighbor_pair_potential(
            md::system const&, md::index, md::index
        )
        {
            return md::softcore_potential<2, 2>{1, 1};
        }
    };

    md::system system;

    system.add_particle().position = {0.5, 0, 0};
    system.add_particle().position = {0, 0.5, 0};
    system.add_particle().position = {0, 0, 0.5};
    std::vector<md::index> reds = {0, 1};

    test_forcefield ff;
    ff.set_subsystem(reds);

    // compute_force does not clear existing force
    std::vector<md::vector> forces = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };
    ff.compute_force(system, forces);

    CHECK(forces[0].x == Approx(1 + 2 * 0.5));
    CHECK(forces[0].y == Approx(2 - 2 * 0.5));
    CHECK(forces[0].z == Approx(3));
    CHECK(forces[1].x == Approx(4 - 2 * 0.5));
    CHECK(forces[1].y == Approx(5 + 2 * 0.5));
    CHECK(forces[1].z == Approx(6));
    CHECK(forces[2].x == Approx(7));
    CHECK(forces[2].y == Approx(8));
    CHECK(forces[2].z == Approx(9));
}

TEST_CASE("make_intra_subsystem_neighbor_pair_forcefield - creates an intra_subsystem_neighbor_pair_forcefield")
{
    md::system system;

    auto ff = md::make_intra_subsystem_neighbor_pair_forcefield(
        md::softcore_potential<2, 2>{1.23, 4.56}
    );
    ff.set_neighbor_distance(4.56);

    auto ndist = ff.intra_subsystem_neighbor_distance(system);
    auto pot = ff.intra_subsystem_neighbor_pair_potential(system, 0, 1);

    using ff_type = decltype(ff);
    using pot_type = decltype(pot);

    CHECK(std::is_base_of<md::intra_subsystem_neighbor_pair_forcefield<ff_type>, ff_type>::value);
    CHECK(std::is_same<pot_type, md::softcore_potential<2, 2>>::value);
    CHECK(pot.energy == 1.23);
    CHECK(pot.diameter == 4.56);
    CHECK(ndist == 4.56);
}
