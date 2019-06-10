#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <md/forcefield/inter_subsystem_pair_forcefield.hpp>

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

TEST_CASE("inter_subsystem_pair_forcefield - computes correct forcefield")
{
    class test_forcefield : public md::inter_subsystem_pair_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential inter_subsystem_pair_potential(
            md::system const&, md::index, md::index
        )
        {
            return md::harmonic_potential{};
        }
    };

    // Particles on two lines.
    md::system system;

    system.add_particle().position = {0, 0, 0};
    system.add_particle().position = {1, 0, 0};
    system.add_particle().position = {2, 0, 0};
    std::vector<md::index> const reds = {0, 1, 2};

    system.add_particle().position = {0, 1, 0};
    system.add_particle().position = {1, 1, 0};
    system.add_particle().position = {2, 1, 0};
    std::vector<md::index> const blacks = {3, 4, 5};

    test_forcefield forcefield;
    forcefield.set_subsystems(reds, blacks);

    // Naive forcefield calculation.
    md::scalar expected_energy = 0;
    std::vector<md::vector> expected_forces(system.particle_count());

    md::array_view<md::point const> positions = system.view_positions();
    md::harmonic_potential potential;

    for (md::index i : reds) {
        for (md::index j : blacks) {
            md::vector const r = positions[i] - positions[j];
            md::vector const force = potential.evaluate_force(r);
            md::scalar const energy = potential.evaluate_energy(r);

            expected_forces[i] += force;
            expected_forces[j] -= force;

            expected_energy += energy;
        }
    }

    // Test.
    md::scalar actual_energy = forcefield.compute_energy(system);
    std::vector<md::vector> actual_forces(system.particle_count());
    forcefield.compute_force(system, actual_forces);

    CHECK(actual_energy == Approx(expected_energy));
    CHECK(max_difference(actual_forces, expected_forces) < 1e-6);
}

TEST_CASE("inter_subsystem_pair_forcefield::compute_force - adds force to array")
{
    class test_forcefield : public md::inter_subsystem_pair_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential inter_subsystem_pair_potential(md::system const&, md::index, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    system.add_particle().position = {1, 0, 0};
    system.add_particle().position = {0, 1, 0};
    std::vector<md::index> const reds = {0};
    std::vector<md::index> const blacks = {1};

    test_forcefield ff;
    ff.set_subsystems(reds, blacks);

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

TEST_CASE("make_inter_subsystem_pair_forcefield - creates an inter_subsystem_pair_forcefield")
{
    md::system system;

    auto ff = md::make_inter_subsystem_pair_forcefield(md::harmonic_potential{1.23});
    auto pot = ff.inter_subsystem_pair_potential(system, 0, 1);

    using ff_type = decltype(ff);
    using pot_type = decltype(pot);

    CHECK(std::is_base_of<md::inter_subsystem_pair_forcefield<ff_type>, ff_type>::value);
    CHECK(std::is_same<pot_type, md::harmonic_potential>::value);
    CHECK(pot.spring_constant == 1.23);
}
