#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/forcefield/plane_surface_forcefield.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <catch.hpp>


TEST_CASE("plane_surface_forcefield - computes inward forcefield")
{
    class inward_forcefield : public md::plane_surface_forcefield<inward_forcefield>
    {
    public:
        md::harmonic_potential plane_inward_potential(md::system const&, md::index) const
        {
            return md::harmonic_potential{};
        }

        md::plane plane_surface(md::system const&) const
        {
            return md::plane{md::vector{1, 0, 0}, md::point{1, 0, 0}};
        }
    };

    md::system system;

    md::scalar const x0 = system.add_particle().position.x = 0.4;
    md::scalar const x1 = system.add_particle().position.x = 0.7;
    md::scalar const x2 = system.add_particle().position.x = 1.0;
    md::scalar const x3 = system.add_particle().position.x = 1.3;
    md::scalar const x4 = system.add_particle().position.x = 1.6;
    (void) x3;
    (void) x4;

    inward_forcefield inward;

    // Energy
    md::scalar const e0 = 0.5 * (x0 - 1) * (x0 - 1);    // inward
    md::scalar const e1 = 0.5 * (x1 - 1) * (x1 - 1);    // inward
    md::scalar const e2 = 0.5 * (x2 - 1) * (x2 - 1);    // on the plane
    md::scalar const e3 = 0;                            // outward
    md::scalar const e4 = 0;                            // outward

    md::scalar const expected_energy = e0 + e1 + e2 + e3 + e4;
    CHECK(inward.compute_energy(system) == Approx(expected_energy));

    // Force
    std::vector<md::vector> forces(system.particle_count());
    inward.compute_force(system, forces);

    CHECK(forces[0].x == Approx(-(x0 - 1)));
    CHECK(forces[1].x == Approx(-(x1 - 1)));
    CHECK(forces[2].x == Approx(-(x2 - 1)));
    CHECK(forces[3].x == Approx(0));
    CHECK(forces[4].x == Approx(0));
}

TEST_CASE("plane_surface_forcefield - computes outward forcefield")
{
    class outward_forcefield : public md::plane_surface_forcefield<outward_forcefield>
    {
    public:
        md::harmonic_potential plane_outward_potential(md::system const&, md::index) const
        {
            return md::harmonic_potential{};
        }

        md::plane plane_surface(md::system const&) const
        {
            return md::plane{md::vector{1, 0, 0}, md::point{1, 0, 0}};
        }
    };

    md::system system;

    md::scalar const x0 = system.add_particle().position.x = 0.4;
    md::scalar const x1 = system.add_particle().position.x = 0.7;
    md::scalar const x2 = system.add_particle().position.x = 1.0;
    md::scalar const x3 = system.add_particle().position.x = 1.3;
    md::scalar const x4 = system.add_particle().position.x = 1.6;
    (void) x0;
    (void) x1;

    outward_forcefield outward;

    // Energy
    md::scalar const e0 = 0;                            // inward
    md::scalar const e1 = 0;                            // inward
    md::scalar const e2 = 0.5 * (x2 - 1) * (x2 - 1);    // on the plane
    md::scalar const e3 = 0.5 * (x3 - 1) * (x3 - 1);    // outward
    md::scalar const e4 = 0.5 * (x4 - 1) * (x4 - 1);    // outward

    md::scalar const expected_energy = e0 + e1 + e2 + e3 + e4;
    CHECK(outward.compute_energy(system) == Approx(expected_energy));

    // Force
    std::vector<md::vector> forces(system.particle_count());
    outward.compute_force(system, forces);

    CHECK(forces[0].x == Approx(0));
    CHECK(forces[1].x == Approx(0));
    CHECK(forces[2].x == Approx(-(x2 - 1)));
    CHECK(forces[3].x == Approx(-(x3 - 1)));
    CHECK(forces[4].x == Approx(-(x4 - 1)));
}

TEST_CASE("plane_surface_forcefield::compute_force - adds force to array")
{
    class outward_forcefield : public md::plane_surface_forcefield<outward_forcefield>
    {
    public:
        md::harmonic_potential plane_outward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;
    outward_forcefield outward;

    // A particle near the xy plane.
    md::scalar const z0 = 0.1;
    md::scalar const f0 = -(z0 - 0);
    system.add_particle().position = {0, 0, z0};

    // compute_force does not clear existing force
    std::vector<md::vector> forces = {
        {1, 2, 3}
    };
    outward.compute_force(system, forces);

    CHECK(forces[0].x == Approx(1));
    CHECK(forces[0].y == Approx(2));
    CHECK(forces[0].z == Approx(3 + f0));
}

TEST_CASE("plane_surface_forcefield::compute_force - collects normal force stats")
{
    class surface_forcefield : public md::plane_surface_forcefield<surface_forcefield>
    {
    public:
        md::harmonic_potential surface_inward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }

        md::harmonic_potential surface_outward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    SECTION("outward force")
    {
        md::system system;
        system.add_particle().position = {0, 0, 1.1};
        system.add_particle().position = {1, 0, 1.2};
        system.add_particle().position = {0, 1, 1.3};

        // xy plane
        surface_forcefield ff;

        std::vector<md::vector> forces(system.particle_count());
        ff.compute_force(system, forces);
        md::scalar const sum_normal = forces[0].z + forces[1].z + forces[2].z;

        CHECK(ff.stats.reaction_force == Approx(-sum_normal));
    }

    SECTION("inward force")
    {
        md::system system;
        system.add_particle().position = {0, 0, 0.9};
        system.add_particle().position = {1, 0, 0.8};
        system.add_particle().position = {0, 1, 0.7};

        // xy plane
        surface_forcefield ff;

        std::vector<md::vector> forces(system.particle_count());
        ff.compute_force(system, forces);
        md::scalar const sum_normal = forces[0].z + forces[1].z + forces[2].z;

        CHECK(ff.stats.reaction_force == Approx(-sum_normal));
    }
}

TEST_CASE("make_plane_inward_forcefield - creates a plane_surface_forcefield")
{
    SECTION("fixed potential")
    {
        md::harmonic_potential potential{1.23};
        auto ff =
            md::make_plane_inward_forcefield(potential)
            .set_plane_surface(
                md::plane{md::vector{0, 0, 1}, md::point{0, 0, 0}}
            );

        md::system system;
        md::harmonic_potential pot = ff.plane_inward_potential(system, 0);
        CHECK(pot.spring_constant == Approx(potential.spring_constant));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::plane_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda potential")
    {
        auto ff =
            md::make_plane_inward_forcefield([](md::index i) {
                return md::harmonic_potential{i * 1.0};
            })
            .set_plane_surface(
                md::plane{md::vector{0, 0, 1}, md::point{0, 0, 0}}
            );

        md::system system;
        md::harmonic_potential pot1 = ff.plane_inward_potential(system, 1);
        md::harmonic_potential pot2 = ff.plane_inward_potential(system, 2);
        CHECK(pot1.spring_constant == Approx(1));
        CHECK(pot2.spring_constant == Approx(2));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::plane_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda plane")
    {
        auto ff =
            md::make_plane_outward_forcefield([](md::index i) {
                return md::harmonic_potential{i * 1.0};
            })
            .set_plane_surface([&] {
                return md::plane{md::vector{1, 2, 3}, md::point{4, 5, 6}};
            });

        md::system system;
        md::plane plane = ff.plane_surface(system);

        CHECK(plane.normal.x == Approx(1));
        CHECK(plane.normal.y == Approx(2));
        CHECK(plane.normal.z == Approx(3));
        CHECK(plane.reference.x == Approx(4));
        CHECK(plane.reference.y == Approx(5));
        CHECK(plane.reference.z == Approx(6));
    }
}

TEST_CASE("make_plane_outward_forcefield - creates a plane_surface_forcefield")
{
    SECTION("fixed potential")
    {
        md::harmonic_potential potential{1.23};
        auto ff =
            md::make_plane_outward_forcefield(potential)
            .set_plane_surface(
                md::plane{md::vector{0, 0, 1}, md::point{0, 0, 0}}
            );

        md::system system;
        md::harmonic_potential pot = ff.plane_outward_potential(system, 0);
        CHECK(pot.spring_constant == Approx(potential.spring_constant));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::plane_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda potential")
    {
        auto ff =
            md::make_plane_outward_forcefield([](md::index i) {
                return md::harmonic_potential{i * 1.0};
            })
            .set_plane_surface(
                md::plane{md::vector{0, 0, 1}, md::point{0, 0, 0}}
            );

        md::system system;
        md::harmonic_potential pot1 = ff.plane_outward_potential(system, 1);
        md::harmonic_potential pot2 = ff.plane_outward_potential(system, 2);
        CHECK(pot1.spring_constant == Approx(1));
        CHECK(pot2.spring_constant == Approx(2));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::plane_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda plane")
    {
        auto ff =
            md::make_plane_outward_forcefield([](md::index i) {
                return md::harmonic_potential{i * 1.0};
            })
            .set_plane_surface([&] {
                return md::plane{md::vector{1, 2, 3}, md::point{4, 5, 6}};
            });

        md::system system;
        md::plane plane = ff.plane_surface(system);

        CHECK(plane.normal.x == Approx(1));
        CHECK(plane.normal.y == Approx(2));
        CHECK(plane.normal.z == Approx(3));
        CHECK(plane.reference.x == Approx(4));
        CHECK(plane.reference.y == Approx(5));
        CHECK(plane.reference.z == Approx(6));
    }
}
