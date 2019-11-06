#include <md/basic_types.hpp>

#include <md/potential/wca_potential.hpp>

#include <catch.hpp>


TEST_CASE("wca_potential - uses defined default parameters")
{
    md::wca_potential pot;

    CHECK(pot.epsilon == 1);
    CHECK(pot.sigma == 1);
}

TEST_CASE("wca_potential - takes zero minimum energy at sigma")
{
    md::wca_potential pot;

    pot.epsilon = 1.23;
    pot.sigma = 4.56;

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

TEST_CASE("wca_potential::evaluate_force - computes correctly directed force")
{
    md::wca_potential pot;

    // Repulsive for r < sigma, zero for r > sigma.
    pot.sigma = 1.5;

    CHECK(pot.evaluate_force({1, 0, 0}).x > 0);
    CHECK(pot.evaluate_force({0, 1, 0}).y > 0);
    CHECK(pot.evaluate_force({0, 0, 1}).z > 0);

    CHECK(pot.evaluate_force({2, 0, 0}).x == Approx(0));
    CHECK(pot.evaluate_force({0, 2, 0}).y == Approx(0));
    CHECK(pot.evaluate_force({0, 0, 2}).z == Approx(0));
}

TEST_CASE("wca_potential - computes correct interaction")
{
    md::wca_potential pot;

    pot.epsilon = 1.2;
    pot.sigma = 3.0;

    md::vector r;

    r = {1, 2, 3};
    CHECK(pot.evaluate_energy(r) == Approx(0));
    CHECK(pot.evaluate_force(r).x == Approx(0));
    CHECK(pot.evaluate_force(r).y == Approx(0));
    CHECK(pot.evaluate_force(r).z == Approx(0));

    r = {0, 1, 2};
    CHECK(pot.evaluate_energy(r) == Approx(28.0179));
    CHECK(pot.evaluate_force(r).x == Approx(0));
    CHECK(pot.evaluate_force(r).y == Approx(81.1590));
    CHECK(pot.evaluate_force(r).z == Approx(162.318));
}
