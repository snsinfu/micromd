#include <md/basic_types.hpp>

#include <md/potential/spring_potential.hpp>

#include <catch.hpp>


TEST_CASE("spring_potential - defaults to centered K=1 spring")
{
    md::spring_potential pot;

    CHECK(pot.spring_constant == 1);
    CHECK(pot.equilibrium_distance == 0);
}

TEST_CASE("spring_potential - is minimum at the equilibrium distance")
{
    md::spring_potential pot;

    md::scalar const k = 1.23;
    md::scalar const b = 4.56;
    pot.spring_constant = k;
    pot.equilibrium_distance = b;

    md::vector const r = md::vector{1, 2, 3}.normalize() * b;
    md::scalar const energy = pot.evaluate_energy(r);
    md::vector const force = pot.evaluate_force(r);

    CHECK(energy == Approx(0).margin(1e-6));
    CHECK(force.x == Approx(0).margin(1e-6));
    CHECK(force.y == Approx(0).margin(1e-6));
    CHECK(force.z == Approx(0).margin(1e-6));
}

TEST_CASE("spring_potential::evaluate_force - computes correctly directed force")
{
    md::spring_potential pot;

    md::scalar const k = 3.21;
    md::scalar const b = 1.23;
    pot.spring_constant = k;
    pot.equilibrium_distance = b;

    // Repulsive for r < b, attractive for r > b.

    CHECK(pot.evaluate_force({1, 0, 0}).x > 0);
    CHECK(pot.evaluate_force({0, 1, 0}).y > 0);
    CHECK(pot.evaluate_force({0, 0, 1}).z > 0);

    CHECK(pot.evaluate_force({2, 0, 0}).x < 0);
    CHECK(pot.evaluate_force({0, 2, 0}).y < 0);
    CHECK(pot.evaluate_force({0, 0, 2}).z < 0);
}

TEST_CASE("spring_potential - computes correct interaction")
{
    md::spring_potential pot;

    md::scalar const k = 1.23;
    md::scalar const b = 4.56;
    pot.spring_constant = k;
    pot.equilibrium_distance = b;

    md::vector const r = {-4, 5, -6};
    md::scalar const energy = pot.evaluate_energy(r);
    md::vector const force = pot.evaluate_force(r);

    CHECK(energy == Approx(k / 2 * (r.norm() - b) * (r.norm() - b)));
    CHECK(force.x == Approx(-k * (r.norm() - b) * r.x / r.norm()));
    CHECK(force.y == Approx(-k * (r.norm() - b) * r.y / r.norm()));
    CHECK(force.z == Approx(-k * (r.norm() - b) * r.z / r.norm()));
}
