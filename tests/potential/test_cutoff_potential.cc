#include <md/basic_types.hpp>
#include <md/potential/cutoff_potential.hpp>
#include <md/potential/lennard_jones_potential.hpp>

#include <catch.hpp>


TEST_CASE("cutoff_potential - cuts potential at configured cutoff distance")
{
    md::lennard_jones_potential pot;
    pot.epsilon = 1;
    pot.sigma = 1;

    md::cutoff_potential<md::lennard_jones_potential> cut = md::apply_cutoff(pot, 2.5);

    SECTION("identity inside cutoff distance")
    {
        md::vector const r1 = {0.8, 0.9, 1.0};

        CHECK(cut.evaluate_energy(r1) == pot.evaluate_energy(r1));
        CHECK(cut.evaluate_force(r1).x == pot.evaluate_force(r1).x);
        CHECK(cut.evaluate_force(r1).y == pot.evaluate_force(r1).y);
        CHECK(cut.evaluate_force(r1).z == pot.evaluate_force(r1).z);
    }

    SECTION("zero at cutoff distance")
    {
        md::vector const r1 = {cut.cutoff_distance, 0, 0};

        CHECK(cut.evaluate_energy(r1) == 0);
        CHECK(cut.evaluate_force(r1).x == 0);
        CHECK(cut.evaluate_force(r1).y == 0);
        CHECK(cut.evaluate_force(r1).z == 0);
    }

    SECTION("zero outside cutoff distance")
    {
        md::vector const r1 = {1.9, 2.0, 2.1};

        CHECK(cut.evaluate_energy(r1) == 0);
        CHECK(cut.evaluate_force(r1).x == 0);
        CHECK(cut.evaluate_force(r1).y == 0);
        CHECK(cut.evaluate_force(r1).z == 0);
    }
}
