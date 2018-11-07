#include <md/basic_types.hpp>

#include <md/potential/softcore_potential.hpp>

#include <catch.hpp>


TEST_CASE("softcore_potential - uses defined default parameters")
{
    md::softcore_potential<3> pot;

    CHECK(pot.overlap_energy == 1);
    CHECK(pot.cutoff_distance == 1);
}

TEST_CASE("softcore_potential - takes maximnum energy at zero")
{
    md::softcore_potential<3> pot;

    pot.overlap_energy = 1.23;
    pot.cutoff_distance = 1;

    md::vector const zero = {};
    CHECK(pot.evaluate_energy(zero) == Approx(1.23));
    CHECK(pot.evaluate_force(zero).x == Approx(0));
    CHECK(pot.evaluate_force(zero).y == Approx(0));
    CHECK(pot.evaluate_force(zero).z == Approx(0));
}

TEST_CASE("softcore_potential - uses specified cutoff distance")
{
    md::softcore_potential<3> pot;

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

TEST_CASE("softcore_potential - supports arbitrary positive exponent")
{
    md::softcore_potential<1> pl1;
    md::softcore_potential<2> pl2;
    md::softcore_potential<3> pl3;
    md::softcore_potential<4> pl4;
    md::softcore_potential<5> pl5;

    md::vector const r = {0.5, 0, 0};

    CHECK(pl1.evaluate_energy(r) == Approx(0.75));
    CHECK(pl2.evaluate_energy(r) == Approx(0.5625));
    CHECK(pl3.evaluate_energy(r) == Approx(0.421875));
    CHECK(pl4.evaluate_energy(r) == Approx(0.31640625));
    CHECK(pl5.evaluate_energy(r) == Approx(0.2373046875));

    CHECK(pl1.evaluate_force(r).x == Approx(1));
    CHECK(pl2.evaluate_force(r).x == Approx(1.5));
    CHECK(pl3.evaluate_force(r).x == Approx(1.6875));
    CHECK(pl4.evaluate_force(r).x == Approx(1.6875));
    CHECK(pl5.evaluate_force(r).x == Approx(1.58203125));
}
