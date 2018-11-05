#include <md/basic_types.hpp>
#include <md/potential/power_law_potential.hpp>

#include <catch.hpp>


TEST_CASE("power_law_potential - uses defined default parameters")
{
    md::power_law_potential<3> pot;

    CHECK(pot.overlap_energy == 1);
    CHECK(pot.cutoff_distance == 1);
}

TEST_CASE("power_law_potential - takes maximnum energy at zero")
{
    md::power_law_potential<3> pot;

    pot.overlap_energy = 1.23;
    pot.cutoff_distance = 1;

    md::vector const zero = {};
    CHECK(pot.evaluate_energy(zero) == Approx(1.23));
    CHECK(pot.evaluate_force(zero).x == Approx(0));
    CHECK(pot.evaluate_force(zero).y == Approx(0));
    CHECK(pot.evaluate_force(zero).z == Approx(0));
}

TEST_CASE("power_law_potential - uses specified cutoff distance")
{
    md::power_law_potential<3> pot;

    pot.overlap_energy = 1.23;
    pot.cutoff_distance = 4.56;

    md::vector const r = {pot.cutoff_distance, 0, 0};

    CHECK(pot.evaluate_energy(r) == Approx(0));
    CHECK(pot.evaluate_force(r).x == Approx(0));
    CHECK(pot.evaluate_force(r).y == Approx(0));
    CHECK(pot.evaluate_force(r).z == Approx(0));

    CHECK(pot.evaluate_energy(1.01 * r) == 0);
    CHECK(pot.evaluate_force(1.01 * r).x == 0);
    CHECK(pot.evaluate_force(1.01 * r).y == 0);
    CHECK(pot.evaluate_force(1.01 * r).z == 0);

    CHECK(pot.evaluate_energy(10 * r) == 0);
    CHECK(pot.evaluate_force(10 * r).x == 0);
    CHECK(pot.evaluate_force(10 * r).y == 0);
    CHECK(pot.evaluate_force(10 * r).z == 0);
}
