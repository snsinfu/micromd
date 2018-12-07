#include <tuple>

#include <md/basic_types.hpp>

#include <md/potential/cosine_bending_potential.hpp>

#include <catch.hpp>


TEST_CASE("cosine_bending_potential - defaults to E=0 cost")
{
    md::cosine_bending_potential pot;

    CHECK(pot.bending_energy == 0);
}

TEST_CASE("cosine_bending_potential - evaluates to zero at zero distance")
{
    md::cosine_bending_potential pot;

    pot.bending_energy = 1.23;

    md::vector const r = {};
    md::vector const s = {};
    md::scalar const energy = pot.evaluate_energy(r, s);
    auto const force = pot.evaluate_force(r, s);

    CHECK(energy == 0);
    CHECK(std::get<0>(force).norm() == 0);
    CHECK(std::get<1>(force).norm() == 0);
    CHECK(std::get<2>(force).norm() == 0);
}

TEST_CASE("cosine_bending_potential - computes correct interaction")
{
    md::cosine_bending_potential pot;

    pot.bending_energy = 1.23;

    md::vector const r = {1, 2, 3};
    md::vector const s = {4, 3, 2};
    md::scalar const energy = pot.evaluate_energy(r, s);
    auto const force = pot.evaluate_force(r, s);

    CHECK(energy == Approx(1.23 * 0.205933));
    CHECK(std::get<0>(force).x == Approx(0.174411));
    CHECK(std::get<0>(force).y == Approx(0.0436028));
    CHECK(std::get<0>(force).z == Approx(-0.0872055));
    CHECK(std::get<1>(force).x == Approx(-0.248085));
    CHECK(std::get<1>(force).y == Approx(-0.0225532));
    CHECK(std::get<1>(force).z == Approx(0.202978));
    CHECK(std::get<2>(force).x == Approx(0.0736736));
    CHECK(std::get<2>(force).y == Approx(-0.0210496));
    CHECK(std::get<2>(force).z == Approx(-0.115773));

}
