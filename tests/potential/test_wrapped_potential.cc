#include <md/basic_types.hpp>
#include <md/potential/wrapped_potential.hpp>
#include <md/potential/scaled_potential.hpp>
#include <md/potential/sum_potential.hpp>

#include <catch.hpp>


TEST_CASE("wrapped_potential - enables arithmetics on custom potential")
{
    struct my_potential
    {
        md::scalar evaluate_energy(md::vector r) const
        {
            return r.squared_norm() / 2;
        }

        md::vector evaluate_force(md::vector r) const
        {
            return r;
        }
    };

    my_potential my1;
    my_potential my2;

    // Check if the wrap_potential result is correctly typed.
    md::wrapped_potential<my_potential> w1 = md::wrap_potential(my1);
    md::wrapped_potential<my_potential> w2 = md::wrap_potential(my2);

    // Allows arithmetics.
    SECTION("energy")
    {
        auto pot = w1 + 0.5 * w2;
        const md::vector r = {1, 2, 3};
        const md::scalar actual = pot.evaluate_energy(r);
        const md::scalar expect = my1.evaluate_energy(r) + 0.5 * my2.evaluate_energy(r);
        CHECK(actual == Approx(expect));
    }

    SECTION("force")
    {
        auto pot = w1 + 0.5 * w2;
        const md::vector r = {1, 2, 3};
        const md::vector actual = pot.evaluate_force(r);
        const md::vector expect = my1.evaluate_force(r) + 0.5 * my2.evaluate_force(r);
        CHECK(actual.x == Approx(expect.x));
        CHECK(actual.y == Approx(expect.y));
        CHECK(actual.z == Approx(expect.z));
    }
}
