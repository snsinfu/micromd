#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <md/forcefield/bonded_pairwise_forcefield.hpp>

#include <catch.hpp>


TEST_CASE("bonded_pairwise_forcefield - computes bond interactions")
{
    class bond_forcefield : public md::bonded_pairwise_forcefield<bond_forcefield>
    {
    public:
        auto bonded_pairwise_potential(md::system const&, md::index, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    md::scalar const x0 = system.add_particle().position.x = 0;
    md::scalar const x1 = system.add_particle().position.x = 1;
    md::scalar const x2 = system.add_particle().position.x = 3;
    md::scalar const x3 = system.add_particle().position.x = 6;
    md::scalar const x4 = system.add_particle().position.x = 10;

    SECTION("no bond")
    {
        bond_forcefield bond;

        // Energy
        CHECK(bond.compute_energy(system) == 0);

        // Force
        std::vector<md::vector> forces(system.particle_count());
        bond.compute_force(system, forces);

        CHECK(forces[0].x == 0);
        CHECK(forces[1].x == 0);
        CHECK(forces[2].x == 0);
        CHECK(forces[3].x == 0);
        CHECK(forces[4].x == 0);
    }

    SECTION("a selected pair")
    {
        // Particle pair 1:2 is bonded.
        md::scalar const expected_energy = 0.5 * (x1 - x2) * (x1 - x2);

        bond_forcefield bond;
        bond.add_bonded_pair(1, 2);

        // Energy
        CHECK(bond.compute_energy(system) == Approx(expected_energy));

        // Force
        std::vector<md::vector> forces(system.particle_count());
        bond.compute_force(system, forces);

        CHECK(forces[0].x == 0);
        CHECK(forces[1].x == Approx(x2 - x1));
        CHECK(forces[2].x == Approx(x1 - x2));
        CHECK(forces[3].x == 0);
        CHECK(forces[4].x == 0);
    }

    SECTION("selected pairs")
    {
        // Particle pairs 0:1, 1:2 and 3:4 are bonded.
        bond_forcefield bond;
        bond.add_bonded_pair(0, 1);
        bond.add_bonded_pair(1, 2);
        bond.add_bonded_pair(3, 4);

        // Energy
        md::scalar const e01 = 0.5 * (x0 - x1) * (x0 - x1);
        md::scalar const e12 = 0.5 * (x1 - x2) * (x1 - x2);
        md::scalar const e34 = 0.5 * (x3 - x4) * (x3 - x4);
        md::scalar const expected_energy = e01 + e12 + e34;

        CHECK(bond.compute_energy(system) == expected_energy);

        // Force
        std::vector<md::vector> forces(system.particle_count());
        bond.compute_force(system, forces);

        CHECK(forces[0].x == Approx(x1 - x0));
        CHECK(forces[1].x == Approx(x0 - x1 + x2 - x1));
        CHECK(forces[2].x == Approx(x1 - x2));
        CHECK(forces[3].x == Approx(x4 - x3));
        CHECK(forces[4].x == Approx(x3 - x4));
    }
}

TEST_CASE("bonded_pairwise_forcefield::add_bonded_pair - returns self")
{
    class test_forcefield : public md::bonded_pairwise_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential bonded_pairwise_potential(md::system const&, md::index, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    test_forcefield test;
    test_forcefield& ref = test.add_bonded_pair(0, 1);

    CHECK(&ref == &test);
}

TEST_CASE("bonded_pairwise_forcefield::compute_force - adds force to array")
{
    class test_forcefield : public md::bonded_pairwise_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential bonded_pairwise_potential(md::system const&, md::index, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    system.add_particle().position = {1, 0, 0};
    system.add_particle().position = {0, 1, 0};

    test_forcefield ff;
    ff.add_bonded_pair(0, 1);

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

TEST_CASE("make_bonded_pairwise_forcefield - creates a bonded_pairwise_forcefield")
{
    md::system system;

    auto ff = md::make_bonded_pairwise_forcefield(md::harmonic_potential{1.23});
    auto pot = ff.bonded_pairwise_potential(system, 0, 1);

    using ff_type = decltype(ff);
    using pot_type = decltype(pot);

    CHECK(std::is_base_of<md::bonded_pairwise_forcefield<ff_type>, ff_type>::value);
    CHECK(std::is_same<pot_type, md::harmonic_potential>::value);
    CHECK(pot.spring_constant == 1.23);
}
