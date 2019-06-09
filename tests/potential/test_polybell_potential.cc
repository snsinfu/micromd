#include <md/basic_types.hpp>

#include <md/potential/polybell_potential.hpp>

#include <catch.hpp>


TEST_CASE("polybell_potential - uses sane default parameters")
{
    md::polybell_potential<4, 4> pot;

    CHECK(pot.overlap_energy == 1);
    CHECK(pot.cutoff_distance == 1);
}

TEST_CASE("polybell_potential - takes maximnum energy at zero")
{
    md::polybell_potential<4, 4> pot;

    pot.overlap_energy = 1.23;
    pot.cutoff_distance = 1;

    md::vector const zero = {};
    CHECK(pot.evaluate_energy(zero) == Approx(1.23));
    CHECK(pot.evaluate_force(zero).x == Approx(0));
    CHECK(pot.evaluate_force(zero).y == Approx(0));
    CHECK(pot.evaluate_force(zero).z == Approx(0));
}

TEST_CASE("polybell_potential - goes to zero at the cutoff distance")
{
    md::polybell_potential<4, 4> pot;

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

TEST_CASE("polybell_potential - computes correct potential")
{
    md::polybell_potential<2, 1> pb1;
    md::polybell_potential<2, 2> pb2;
    md::polybell_potential<3, 3> pb3;
    md::polybell_potential<4, 4> pb4;
    md::polybell_potential<5, 5> pb5;

    md::vector const r = {0.5, 0, 0};

    CHECK(pb1.evaluate_energy(r) == Approx(0.75));
    CHECK(pb2.evaluate_energy(r) == Approx(0.5625));
    CHECK(pb3.evaluate_energy(r) == Approx(0.669921875));
    CHECK(pb4.evaluate_energy(r) == Approx(0.7724761962890625));
    CHECK(pb5.evaluate_energy(r) == Approx(0.8532151877880096));

    CHECK(pb1.evaluate_force(r).x == Approx(1.0));
    CHECK(pb2.evaluate_force(r).x == Approx(1.5));
    CHECK(pb3.evaluate_force(r).x == Approx(1.72265625));
    CHECK(pb4.evaluate_force(r).x == Approx(1.64794921875));
    CHECK(pb5.evaluate_force(r).x == Approx(1.3761535286903381));
}
