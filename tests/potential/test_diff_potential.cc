#include <md/basic_types.hpp>
#include <md/potential/harmonic_potential.hpp>
#include <md/potential/lennard_jones_potential.hpp>
#include <md/potential/spring_potential.hpp>
#include <md/potential/diff_potential.hpp>

#include <catch.hpp>


TEST_CASE("diff_potential - computes the diff potential")
{
    md::harmonic_potential harmonic{0.1};
    md::spring_potential spring{0.2, 0.3};

    using diff_type = md::diff_potential<md::harmonic_potential, md::spring_potential>;
    const diff_type diff = harmonic - spring;

    SECTION("energy")
    {
        const md::vector r = {1, 2, 3};
        const md::scalar actual = diff.evaluate_energy(r);
        const md::scalar expect = harmonic.evaluate_energy(r) - spring.evaluate_energy(r);
        CHECK(actual == Approx(expect));
    }

    SECTION("force")
    {
        const md::vector r = {1, 2, 3};
        const md::vector actual = diff.evaluate_force(r);
        const md::vector expect = harmonic.evaluate_force(r) - spring.evaluate_force(r);
        CHECK(actual.x == Approx(expect.x));
        CHECK(actual.y == Approx(expect.y));
        CHECK(actual.z == Approx(expect.z));
    }
}

TEST_CASE("diff_potential - can be an operand of another diff")
{
    md::harmonic_potential harmonic{0.1};
    md::spring_potential spring{0.2, 0.3};
    md::lennard_jones_potential lj{0.4, 0.5};

    using diff_type = md::diff_potential<
        md::diff_potential<md::harmonic_potential, md::spring_potential>,
        md::lennard_jones_potential
    >;
    const diff_type diff = harmonic - spring - lj;

    SECTION("energy")
    {
        const md::vector r = {1, 2, 3};
        const md::scalar actual = diff.evaluate_energy(r);
        const md::scalar expect =
            harmonic.evaluate_energy(r) - spring.evaluate_energy(r) - lj.evaluate_energy(r);
        CHECK(actual == Approx(expect));
    }

    SECTION("force")
    {
        const md::vector r = {1, 2, 3};
        const md::vector actual = diff.evaluate_force(r);
        const md::vector expect =
            harmonic.evaluate_force(r) - spring.evaluate_force(r) - lj.evaluate_force(r);
        CHECK(actual.x == Approx(expect.x));
        CHECK(actual.y == Approx(expect.y));
        CHECK(actual.z == Approx(expect.z));
    }
}
