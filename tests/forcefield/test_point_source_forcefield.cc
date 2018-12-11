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
    md::system system;

    auto ff = md::make_point_source_forcefield(md::harmonic_potential{1.23});
    ff.set_point_source(md::point{1, 2, 3});

    auto pot = ff.point_source_potential(system, 0);

    using ff_type = decltype(ff);
    using pot_type = decltype(pot);

    CHECK(std::is_base_of<md::point_source_forcefield<ff_type>, ff_type>::value);
    CHECK(std::is_same<pot_type, md::harmonic_potential>::value);
    CHECK(pot.spring_constant == 1.23);
}
