#include <md/basic_types.hpp>

#include <md/potential/soft_lennard_jones_potential.hpp>

#include <catch.hpp>


TEST_CASE("soft_lennard_jones_potential - uses sane default paramters")
{
    md::soft_lennard_jones_potential pot;

    CHECK(pot.epsilon == 1);
    CHECK(pot.sigma == 1);
    CHECK(pot.softness == 0.1);
}

TEST_CASE("soft_lennard_jones_potential - takes minimum energy at sigma")
{
    md::soft_lennard_jones_potential pot;

    pot.epsilon = 1.23;
    pot.sigma = 4.56;
    pot.softness = 7.89;

    SECTION("x direction")
    {
        md::vector const r = {pot.sigma, 0, 0};
        md::scalar const energy = pot.evaluate_energy(r);
        md::vector const force = pot.evaluate_force(r);

        CHECK(energy == Approx(0));
        CHECK(force.x == Approx(0));
        CHECK(force.y == Approx(0));
        CHECK(force.z == Approx(0));
    }

    SECTION("y direction")
    {
        md::vector const r = {0, pot.sigma, 0};
        md::scalar const energy = pot.evaluate_energy(r);
        md::vector const force = pot.evaluate_force(r);

        CHECK(energy == Approx(0));
        CHECK(force.x == Approx(0));
        CHECK(force.y == Approx(0));
        CHECK(force.z == Approx(0));
    }

    SECTION("z direction")
    {
        md::vector const r = {0, 0, pot.sigma};
        md::scalar const energy = pot.evaluate_energy(r);
        md::vector const force = pot.evaluate_force(r);

        CHECK(energy == Approx(0));
        CHECK(force.x == Approx(0));
        CHECK(force.y == Approx(0));
        CHECK(force.z == Approx(0));
    }
}

TEST_CASE("soft_lennard_jones_potential - has a finite extremum at zero distance")
{
    md::soft_lennard_jones_potential pot;

    pot.epsilon = 1.23;
    pot.sigma = 4.56;
    pot.softness = 7.89;

    md::vector const r = {0, 0, 0};
    md::scalar const energy = pot.evaluate_energy(r);
    md::vector const force = pot.evaluate_force(r);

    CHECK(energy == Approx(pot.epsilon / (pot.softness * pot.softness)));
    CHECK(force.x == Approx(0));
    CHECK(force.y == Approx(0));
    CHECK(force.z == Approx(0));
}

TEST_CASE("soft_lennard_jones_potential::evaluate_force - computes correctly directed force")
{
    md::soft_lennard_jones_potential pot;

    // Repulsive for r < sigma, attractive for r > sigma.
    pot.sigma = 1.5;

    CHECK(pot.evaluate_force({1, 0, 0}).x > 0);
    CHECK(pot.evaluate_force({0, 1, 0}).y > 0);
    CHECK(pot.evaluate_force({0, 0, 1}).z > 0);

    CHECK(pot.evaluate_force({2, 0, 0}).x < 0);
    CHECK(pot.evaluate_force({0, 2, 0}).y < 0);
    CHECK(pot.evaluate_force({0, 0, 2}).z < 0);
}

TEST_CASE("soft_lennard_jones_potential - computes correct interaction")
{
    md::soft_lennard_jones_potential pot;

    pot.epsilon = 1.2;
    pot.sigma = 3.0;
    pot.softness = 0.1;

    md::vector r;

    r = {1, 2, 3};
    CHECK(pot.evaluate_energy(r) == Approx(0.614028));
    CHECK(pot.evaluate_force(r).x == Approx(-0.788394));
    CHECK(pot.evaluate_force(r).y == Approx(-1.576788));
    CHECK(pot.evaluate_force(r).z == Approx(-2.365183));

    r = {0, 1, 2};
    CHECK(pot.evaluate_energy(r) == Approx(11.177985));
    CHECK(pot.evaluate_force(r).x == Approx(0));
    CHECK(pot.evaluate_force(r).y == Approx(6.107177));
    CHECK(pot.evaluate_force(r).z == Approx(12.214355));
}
