#include <md/basic_types.hpp>

#include <md/potential/constant_potential.hpp>

#include <catch.hpp>


TEST_CASE("constant_potential - defaults to zero")
{
    md::constant_potential pot;

    CHECK(pot.energy == 0);
}

TEST_CASE("constant_potential::evaluate_energy - returns configured energy")
{
    md::constant_potential pot;

    pot.energy = 100;

    CHECK(pot.evaluate_energy({0, 0, 0}) == 100);
    CHECK(pot.evaluate_energy({1, 0, 0}) == 100);
    CHECK(pot.evaluate_energy({9, 0, 0}) == 100);
}

TEST_CASE("constant_potential::evaluate_force - returns zero")
{
    md::constant_potential pot;

    pot.energy = 100;

    CHECK(pot.evaluate_force({0, 0, 0}).norm() == 0);
    CHECK(pot.evaluate_force({1, 0, 0}).norm() == 0);
    CHECK(pot.evaluate_force({9, 0, 0}).norm() == 0);
}
