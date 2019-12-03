#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <md/forcefield/bruteforce_pairwise_forcefield.hpp>

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

TEST_CASE("bruteforce_pairwise_forcefield - computes correct forcefield")
{
    class test_forcefield : public md::bruteforce_pairwise_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential bruteforce_pairwise_potential(md::system const&, md::index, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    // Particles on a 5x5x5 grid
    md::system system;

    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            for (int z = -2; z <= 2; z++) {
                system.add_particle().position = {
                    x * 0.5,
                    y * 0.5,
                    z * 0.5
                };
            }
        }
    }

    test_forcefield forcefield;
    md::harmonic_potential potential;

    std::vector<md::vector> actual_forces(system.particle_count());
    std::vector<md::vector> expected_forces(system.particle_count());

    forcefield.compute_force(system, actual_forces);

    md::scalar actual_energy = forcefield.compute_energy(system);
    md::scalar expected_energy = 0;

    md::array_view<md::point const> positions = system.view_positions();

    for (md::index i = 0; i < positions.size(); i++) {
        for (md::index j = 0; j < i; j++) {
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

TEST_CASE("bruteforce_pairwise_forcefield::compute_force - adds force to array")
{
    class test_forcefield : public md::bruteforce_pairwise_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential bruteforce_pairwise_potential(md::system const&, md::index, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    system.add_particle().position = {1, 0, 0};
    system.add_particle().position = {0, 1, 0};

    test_forcefield ff;

    // compute_force does not clear existing force
    std::vector<md::vector> forces = {
        {1, 2, 3},
        {4, 5, 6}
    };
    ff.compute_force(system, forces);

    CHECK(forces[0].x == Approx(1 - 1));
    CHECK(forces[0].y == Approx(2 + 1));
    CHECK(forces[0].z == Approx(3));
    CHECK(forces[1].x == Approx(4 + 1));
    CHECK(forces[1].y == Approx(5 - 1));
    CHECK(forces[1].z == Approx(6));
}

TEST_CASE("make_bruteforce_pairwise_forcefield - creates an bruteforce_pairwise_forcefield")
{
    md::system system;

    auto ff = md::make_bruteforce_pairwise_forcefield(md::harmonic_potential{1.23});
    auto pot = ff.bruteforce_pairwise_potential(system, 0, 1);

    using ff_type = decltype(ff);
    using pot_type = decltype(pot);

    CHECK(std::is_base_of<md::bruteforce_pairwise_forcefield<ff_type>, ff_type>::value);
    CHECK(std::is_same<pot_type, md::harmonic_potential>::value);
    CHECK(pot.spring_constant == 1.23);
}
