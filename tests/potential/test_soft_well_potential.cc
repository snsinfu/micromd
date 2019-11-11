#include <md/basic_types.hpp>
#include <md/potential/soft_well_potential.hpp>

#include <catch.hpp>


TEST_CASE("soft_well_potential - uses sane default paramters")
{
    md::soft_well_potential<2> pot;

    CHECK(pot.energy == 1);
    CHECK(pot.decay_distance == 1);
}

TEST_CASE("soft_well_potential - is minimum at zero distance")
{
    md::soft_well_potential<2> pot;
    pot.energy = 1.23;
    pot.decay_distance = 4.56;

    md::vector const r = {0, 0, 0};
    md::scalar const energy = pot.evaluate_energy(r);
    md::vector const force = pot.evaluate_force(r);

    CHECK(energy == Approx(-pot.energy));
    CHECK(force.norm() == Approx(0));
}

TEST_CASE("soft_well_potential - decays to zero")
{
    md::soft_well_potential<2> pot;
    pot.energy = 1.23;
    pot.decay_distance = 4.56;

    md::vector const r = {100, 0, 0};
    md::scalar const energy = pot.evaluate_energy(r);
    md::vector const force = pot.evaluate_force(r);

    CHECK(energy < 0);
    CHECK(energy == Approx(0).margin(0.01));

    CHECK(force.dot(r) < 0); // Attractive
    CHECK(force.norm() == Approx(0).margin(0.01));
}

TEST_CASE("soft_well_potential - halves at decay_distance")
{
    md::soft_well_potential<2> pot;
    pot.energy = 1.23;
    pot.decay_distance = 4.56;

    md::scalar const energy = pot.evaluate_energy({pot.decay_distance, 0, 0});
    CHECK(energy == Approx(-pot.energy / 2));
}

TEST_CASE("soft_well_potential - computes correct values")
{
    SECTION("2nd order")
    {
        md::soft_well_potential<2> pot;
        pot.energy = 1.23;
        pot.decay_distance = 4.56;

        md::vector const r1 = {1.2, 3.4, 5.6};
        md::scalar const e1 = pot.evaluate_energy(r1);
        md::vector const f1 = pot.evaluate_force(r1);

        CHECK(e1 == Approx(-0.39255126347584784));
        CHECK(f1.x == Approx(-0.01446003033358149));
        CHECK(f1.y == Approx(-0.040970085945147554));
        CHECK(f1.z == Approx(-0.06748014155671363));
    }

    SECTION("6th order")
    {
        md::soft_well_potential<6> pot;
        pot.energy = 1.23;
        pot.decay_distance = 4.56;

        md::vector const r1 = {1.2, 3.4, 5.6};
        md::scalar const e1 = pot.evaluate_energy(r1);
        md::vector const f1 = pot.evaluate_force(r1);

        CHECK(e1 == Approx(-0.11485401188866758));
        CHECK(f1.x == Approx(-0.016901052221443922));
        CHECK(f1.y == Approx(-0.047886314627424445));
        CHECK(f1.z == Approx(-0.07887157703340497));
    }
}
