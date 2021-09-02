#include <tuple>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/forcefield/bonded_triplewise_forcefield.hpp>

#include <catch.hpp>


namespace
{
    // u(r,s) = dot(r, s)
    struct dot_potential
    {
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

TEST_CASE("bonded_triplewise_forcefield - computes interactions among triplets")
{
    class chain_forcefield : public md::bonded_triplewise_forcefield<chain_forcefield>
    {
    public:
        auto bonded_triplewise_potential(
            md::system const&, md::index, md::index, md::index
        ) const
        {
            return dot_potential{};
        }
    };

    md::system system;

    md::point const p0 = system.add_particle().position = {0, 0, 0};
    md::point const p1 = system.add_particle().position = {1, 0, 0};
    md::point const p2 = system.add_particle().position = {3, 0, 0};
    md::point const p3 = system.add_particle().position = {6, 0, 0};
    md::point const p4 = system.add_particle().position = {10, 0, 0};

    SECTION("no triples")
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

    SECTION("selected triples")
    {
        chain_forcefield chain;

        // 0:1:2 and 2:3:4
        chain.add_bonded_triple(0, 1, 2);
        chain.add_bonded_triple(2, 3, 4);

        md::vector const r01 = p0 - p1;
        md::vector const r12 = p1 - p2;
        md::vector const r23 = p2 - p3;
        md::vector const r34 = p3 - p4;

        dot_potential potential;

        // Energy
        md::scalar const expected_energy =
            potential.evaluate_energy(r01, r12) +
            potential.evaluate_energy(r23, r34);
        CHECK(chain.compute_energy(system) == Approx(expected_energy));

        // Force
        std::vector<md::vector> forces(system.particle_count());
        chain.compute_force(system, forces);

        auto const f012 = potential.evaluate_force(r01, r12);
        auto const f234 = potential.evaluate_force(r23, r34);
        md::vector const f0 = std::get<0>(f012);
        md::vector const f1 = std::get<1>(f012);
        md::vector const f2 = std::get<2>(f012) + std::get<0>(f234);
        md::vector const f3 =                     std::get<1>(f234);
        md::vector const f4 =                     std::get<2>(f234);

        CHECK(forces[0].x == Approx(f0.x));
        CHECK(forces[1].x == Approx(f1.x));
        CHECK(forces[2].x == Approx(f2.x));
        CHECK(forces[3].x == Approx(f3.x));
        CHECK(forces[4].x == Approx(f4.x));
    }

    SECTION("range of triples")
    {
        chain_forcefield chain;

        // 1:2:3 and 2:3:4
        chain.add_bonded_range(1, 5);

        md::vector const r12 = p1 - p2;
        md::vector const r23 = p2 - p3;
        md::vector const r34 = p3 - p4;

        dot_potential potential;

        // Energy
        md::scalar const expected_energy =
            potential.evaluate_energy(r12, r23) +
            potential.evaluate_energy(r23, r34);
        CHECK(chain.compute_energy(system) == Approx(expected_energy));

        // Force
        std::vector<md::vector> forces(system.particle_count());
        chain.compute_force(system, forces);

        auto const f123 = potential.evaluate_force(r12, r23);
        auto const f234 = potential.evaluate_force(r23, r34);
        md::vector const f0 = {};
        md::vector const f1 = std::get<0>(f123);
        md::vector const f2 = std::get<1>(f123) + std::get<0>(f234);
        md::vector const f3 = std::get<2>(f123) + std::get<1>(f234);
        md::vector const f4 =                     std::get<2>(f234);

        CHECK(forces[0].x == Approx(f0.x));
        CHECK(forces[1].x == Approx(f1.x));
        CHECK(forces[2].x == Approx(f2.x));
        CHECK(forces[3].x == Approx(f3.x));
        CHECK(forces[4].x == Approx(f4.x));
    }
}

TEST_CASE("bonded_triplewise_forcefield::add_bonded_triple - returns self")
{
    class test_forcefield : public md::bonded_triplewise_forcefield<test_forcefield>
    {
    public:
        dot_potential bonded_triplewise_potential(
            md::system const&, md::index, md::index, md::index
        ) const
        {
            return dot_potential{};
        }
    };

    test_forcefield test;
    test_forcefield& ref = test.add_bonded_triple(0, 1, 2);

    CHECK(&ref == &test);
}

TEST_CASE("bonded_triplewise_forcefield::add_bonded_range - returns self")
{
    class test_forcefield : public md::bonded_triplewise_forcefield<test_forcefield>
    {
    public:
        dot_potential bonded_triplewise_potential(
            md::system const&, md::index, md::index, md::index
        ) const
        {
            return dot_potential{};
        }
    };

    test_forcefield test;
    test_forcefield& ref = test.add_bonded_range(0, 10);

    CHECK(&ref == &test);
}

TEST_CASE("bonded_triplewise_forcefield::compute_force - adds force to array")
{
    class test_forcefield : public md::bonded_triplewise_forcefield<test_forcefield>
    {
    public:
        dot_potential bonded_triplewise_potential(
            md::system const&, md::index, md::index, md::index
        ) const
        {
            return dot_potential{};
        }
    };

    md::system system;

    md::point const p0 = system.add_particle().position = {1, 0, 0};
    md::point const p1 = system.add_particle().position = {0, 1, 0};
    md::point const p2 = system.add_particle().position = {0, 0, 1};

    test_forcefield ff;
    ff.add_bonded_triple(0, 1, 2);

    // compute_force does not clear existing force
    std::vector<md::vector> forces = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };
    ff.compute_force(system, forces);

    CHECK(forces[0].x == Approx(1 - (p1 - p2).x));
    CHECK(forces[0].y == Approx(2 - (p1 - p2).y));
    CHECK(forces[0].z == Approx(3 - (p1 - p2).z));

    CHECK(forces[1].x == Approx(4 + (p1 - p2).x - (p0 - p1).x));
    CHECK(forces[1].y == Approx(5 + (p1 - p2).y - (p0 - p1).y));
    CHECK(forces[1].z == Approx(6 + (p1 - p2).z - (p0 - p1).z));

    CHECK(forces[2].x == Approx(7 + (p0 - p1).x));
    CHECK(forces[2].y == Approx(8 + (p0 - p1).y));
    CHECK(forces[2].z == Approx(9 + (p0 - p1).z));
}

TEST_CASE("make_bonded_triplewise_forcefield - creates a bonded_triplewise_forcefield")
{
    SECTION("potential")
    {
        auto ff = md::make_bonded_triplewise_forcefield(dot_potential{});

        md::system system;
        ff.bonded_triplewise_potential(system, 0, 1, 2);

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::bonded_triplewise_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda (i, j, k)")
    {
        auto ff = md::make_bonded_triplewise_forcefield(
            [](md::index, md::index, md::index) {
                return dot_potential{};
            }
        );

        md::system system;
        ff.bonded_triplewise_potential(system, 0, 1, 2);
    }

    SECTION("lambda (system, i, j, k)")
    {
        auto ff = md::make_bonded_triplewise_forcefield(
            [](md::system const&, md::index, md::index, md::index) {
                return dot_potential{};
            }
        );

        md::system system;
        ff.bonded_triplewise_potential(system, 0, 1, 2);
    }
}
