#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <md/forcefield/point_source_forcefield.hpp>

#include <catch.hpp>


TEST_CASE("point_source_forcefield - computes correct forcefield")
{
    class test_forcefield : public md::point_source_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential point_source_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;
    test_forcefield forcefield;

    md::point const p0 = system.add_particle().position = {0, 1, 2};
    md::point const p1 = system.add_particle().position = {3, 4, 5};
    md::point const p2 = system.add_particle().position = {6, 7, 8};
    md::point const ps = {9, 0, 1};
    forcefield.set_point_source(ps);

    md::scalar const expected_energy =
        0.5 * (p0 - ps).squared_norm() +
        0.5 * (p1 - ps).squared_norm() +
        0.5 * (p2 - ps).squared_norm();

    CHECK(forcefield.compute_energy(system) == Approx(expected_energy));

    // Force
    std::vector<md::vector> forces(system.particle_count());
    forcefield.compute_force(system, forces);

    CHECK((forces[0] - (ps - p0)).norm() == Approx(0).margin(1e-6));
    CHECK((forces[1] - (ps - p1)).norm() == Approx(0).margin(1e-6));
    CHECK((forces[2] - (ps - p2)).norm() == Approx(0).margin(1e-6));
}

TEST_CASE("point_source_forcefield::set_point_source - returns self")
{
    class test_forcefield : public md::point_source_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential point_source_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    test_forcefield test;
    test_forcefield& ref = test.set_point_source(md::point{0, 0, 0});

    CHECK(&ref == &test);
}

TEST_CASE("point_source_forcefield::set_point_source_targets - returns self")
{
    class test_forcefield : public md::point_source_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential point_source_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    std::vector<md::index> targets;
    test_forcefield test;
    test_forcefield& ref = test.set_point_source_targets(targets);

    CHECK(&ref == &test);
}

TEST_CASE("point_source_forcefield::compute_force - adds force to array")
{
    class test_forcefield : public md::point_source_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential point_source_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    md::point const p0 = system.add_particle().position = {1.1, 0, 0};

    test_forcefield test;
    test.set_point_source(md::point{0, 0, 0});

    // compute_force does not clear existing force
    std::vector<md::vector> forces = {
        {1, 2, 3}
    };
    test.compute_force(system, forces);

    CHECK(forces[0].x == Approx(1 - p0.x));
    CHECK(forces[0].y == Approx(2 - p0.y));
    CHECK(forces[0].z == Approx(3 - p0.z));
}

TEST_CASE("make_point_source_forcefield - creates a point_source_forcefield")
{
    SECTION("fixed potential")
    {
        auto ff = md::make_point_source_forcefield(md::harmonic_potential{1.23});
        ff.set_point_source(md::point{1, 2, 3});

        md::system system;
        md::harmonic_potential pot = ff.point_source_potential(system, 0);
        CHECK(pot.spring_constant == 1.23);

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::point_source_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda potential")
    {
        auto ff = md::make_point_source_forcefield([](md::index i) {
            return md::harmonic_potential{i * 1.0};
        });
        ff.set_point_source(md::point{1, 2, 3});

        md::system system;
        md::harmonic_potential pot1 = ff.point_source_potential(system, 1);
        md::harmonic_potential pot2 = ff.point_source_potential(system, 2);
        CHECK(pot1.spring_constant == Approx(1));
        CHECK(pot2.spring_constant == Approx(2));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::point_source_forcefield<ff_type>, ff_type>::value);
    }
}

TEST_CASE("point_source_forcefield::set_point_source_targets - selects forcefield targets")
{
    md::system system;
    system.add_particle().position = {1, 2, 3};
    system.add_particle().position = {4, 5, 6};
    system.add_particle().position = {7, 8, 9};
    system.add_particle().position = {0, 1, 2};
    system.add_particle().position = {3, 4, 5};

    md::point const source = {0, 0, 0};
    std::vector<md::index> const targets = {1, 3, 4};

    md::harmonic_potential potential;
    system.add_forcefield(
        md::make_point_source_forcefield(potential)
        .set_point_source(source)
        .set_point_source_targets(targets)
    );

    auto const positions = system.view_positions();

    SECTION("energy")
    {
        md::scalar expected_energy = 0;
        for (md::index const i : targets) {
            expected_energy += potential.evaluate_energy(positions[i] - source);
        }

        md::scalar const actual_energy = system.compute_energy();

        CHECK(actual_energy == Approx(expected_energy));
    }

    SECTION("force")
    {
        std::vector<md::vector> expected_forces(system.particle_count());
        for (md::index const i : targets) {
            expected_forces[i] = potential.evaluate_force(positions[i] - source);
        }

        std::vector<md::vector> actual_forces(system.particle_count());
        system.compute_force(actual_forces);

        CHECK((actual_forces[0] - expected_forces[0]).norm() == Approx(0));
        CHECK((actual_forces[1] - expected_forces[1]).norm() == Approx(0));
        CHECK((actual_forces[2] - expected_forces[2]).norm() == Approx(0));
        CHECK((actual_forces[3] - expected_forces[3]).norm() == Approx(0));
        CHECK((actual_forces[4] - expected_forces[4]).norm() == Approx(0));
    }
}
