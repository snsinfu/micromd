#include <md/basic_types.hpp>
#include <md/potential/harmonic_potential.hpp>
#include <md/potential/lennard_jones_potential.hpp>
#include <md/potential/spring_potential.hpp>
#include <md/potential/scaled_potential.hpp>

#include <catch.hpp>


TEST_CASE("scaled_potential - computes scalar-multiplied potential")
{
    md::harmonic_potential harmonic{0.1};

    using scaled_type = md::scaled_potential<md::harmonic_potential>;
    const md::scalar factor = 5;
    const scaled_type scaled = factor * harmonic;

    SECTION("energy")
    {
        const md::vector r = {1, 2, 3};
        const md::scalar actual = scaled.evaluate_energy(r);
        const md::scalar expect = factor * harmonic.evaluate_energy(r);
        CHECK(actual == Approx(expect));
    }

    SECTION("force")
    {
        const md::vector r = {1, 2, 3};
        const md::vector actual = scaled.evaluate_force(r);
        const md::vector expect = factor * harmonic.evaluate_force(r);
        CHECK(actual.x == Approx(expect.x));
        CHECK(actual.y == Approx(expect.y));
        CHECK(actual.z == Approx(expect.z));
    }
}
