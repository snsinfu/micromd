#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <md/forcefield/sequential_pair_forcefield.hpp>

#include <catch.hpp>


TEST_CASE("sequential_pair_forcefield - computes bond interactions")
{
    class bond_forcefield : public md::sequential_pair_forcefield<bond_forcefield>
    {
    public:
        auto sequential_pair_potential(md::system const&, md::index, md::index)
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

    SECTION("a bonded segment")
    {
        // Particle pairs 1:2 and 2:3 are bonded.
        md::scalar const e12 = 0.5 * (x1 - x2) * (x1 - x2);
        md::scalar const e23 = 0.5 * (x2 - x3) * (x2 - x3);
        md::scalar const expected_energy = e12 + e23;

        bond_forcefield bond;
        bond.add_segment(1, 3);

        // Energy
        CHECK(bond.compute_energy(system) == Approx(expected_energy));

        // Force
        std::vector<md::vector> forces(system.particle_count());
        bond.compute_force(system, forces);

        CHECK(forces[0].x == 0);
        CHECK(forces[1].x == Approx(x2 - x1));
        CHECK(forces[2].x == Approx(x1 - x2 + x3 - x2));
        CHECK(forces[3].x == Approx(x2 - x3));
        CHECK(forces[4].x == 0);
    }

    SECTION("two bonded segments")
    {
        // Particle pairs 0:1, 1:2 and 3:4 are bonded.
        bond_forcefield bond;
        bond.add_segment(0, 2);
        bond.add_segment(3, 4);

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
