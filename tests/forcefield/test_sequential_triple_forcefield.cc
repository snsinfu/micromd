#include <tuple>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/forcefield/sequential_triple_forcefield.hpp>

#include <catch.hpp>


namespace
{
    // u(r,s) = K dot(r, s)
    struct dot_potential
    {
        md::scalar K;

        md::scalar evaluate_energy(md::vector r, md::vector s) const
        {
            return r.dot(s);
        }

        auto evaluate_force(md::vector r, md::vector s) const
        {
            return std::make_tuple(-s, s - r, r);
        }
    };
}

TEST_CASE("sequential_triple_forcefield - computes interactions within triplets")
{
    class chain_forcefield : public md::sequential_triple_forcefield<chain_forcefield>
    {
    public:
        auto sequential_triple_potential(md::system const&, md::index, md::index, md::index)
        {
            return dot_potential{};
        }
    };

    md::system system;

    md::scalar const x0 = system.add_particle().position.x = 0;
    md::scalar const x1 = system.add_particle().position.x = 1;
    md::scalar const x2 = system.add_particle().position.x = 3;
    md::scalar const x3 = system.add_particle().position.x = 6;
    md::scalar const x4 = system.add_particle().position.x = 10;

    SECTION("no segment")
    {
        chain_forcefield chain;

        // Energy
        CHECK(chain.compute_energy(system) == 0);

        // Force
        std::vector<md::vector> forces(system.particle_count());
        chain.compute_force(system, forces);

        CHECK(forces[0].x == 0);
        CHECK(forces[1].x == 0);
        CHECK(forces[2].x == 0);
        CHECK(forces[3].x == 0);
        CHECK(forces[4].x == 0);
    }

    SECTION("one segment")
    {
        // Particle triples 0:1:2 and 1:2:3 are considered.
        md::scalar const e012 = (x0 - x1) * (x1 - x2);
        md::scalar const e123 = (x1 - x2) * (x2 - x3);
        md::scalar const expected_energy = e012 + e123;

        chain_forcefield chain;
        chain.add_segment(0, 3);

        // Energy
        CHECK(chain.compute_energy(system) == Approx(expected_energy));

        // Force
        std::vector<md::vector> forces(system.particle_count());
        chain.compute_force(system, forces);

        CHECK(forces[0].x == Approx(-x1 + x2));
        CHECK(forces[1].x == Approx(-x0 + 2 * x1 - 2 * x2 + x3));
        CHECK(forces[2].x == Approx(x0 - 2 * x1 + 2 * x2 - x3));
        CHECK(forces[3].x == Approx(x1 - x2));
        CHECK(forces[4].x == 0);
    }
}

TEST_CASE("sequential_triple_forcefield::add_segment - returns self")
{
    class test_forcefield : public md::sequential_triple_forcefield<test_forcefield>
    {
    public:
        dot_potential sequential_triple_potential(
            md::system const&, md::index, md::index, md::index
        ) const
        {
            return dot_potential{};
        }
    };

    test_forcefield test;
    test_forcefield& ref = test.add_segment(0, 9);

    CHECK(&ref == &test);
}

TEST_CASE("make_sequential_triple_forcefield - creates a sequential_triple_forcefield")
{
    md::system system;

    auto ff = md::make_sequential_triple_forcefield(dot_potential{1.23});
    auto pot = ff.sequential_triple_potential(system, 0, 1, 2);

    using ff_type = decltype(ff);
    using pot_type = decltype(pot);

    CHECK(std::is_base_of<md::sequential_triple_forcefield<ff_type>, ff_type>::value);
    CHECK(std::is_same<pot_type, dot_potential>::value);
    CHECK(pot.K == 1.23);
}

