#include <md/basic_types.hpp>

#include <md/potential/harmonic_potential.hpp>

#include <catch.hpp>


TEST_CASE("harmonic_potential - defaults to K=1 spring")
{
    md::harmonic_potential pot;

    CHECK(pot.spring_constant == 1);
}

TEST_CASE("harmonic_potential - evaluates to zero at zero distance")
{
    md::harmonic_potential pot;

    md::vector const r = {};
    md::scalar const energy = pot.evaluate_energy(r);
    md::vector const force = pot.evaluate_force(r);

    CHECK(energy == 0);
    CHECK(force.x == 0);
    CHECK(force.y == 0);
    CHECK(force.z == 0);
}

TEST_CASE("harmonic_potential - computes correct interaction")
{
    md::harmonic_potential pot;

    pot.spring_constant = 1.23;

    md::vector const r = {-4, 5, -6};
    md::scalar const energy = pot.evaluate_energy(r);
    md::vector const force = pot.evaluate_force(r);

    CHECK(energy == Approx(0.5 * pot.spring_constant * r.squared_norm()));
    CHECK(force.x == Approx(-pot.spring_constant * r.x));
    CHECK(force.y == Approx(-pot.spring_constant * r.y));
    CHECK(force.z == Approx(-pot.spring_constant * r.z));
}
