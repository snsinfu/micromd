#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <md/forcefield/sphere_surface_forcefield.hpp>

#include <catch.hpp>


TEST_CASE("sphere::implicit - returns negative inside, positive outside")
{
    md::sphere sphere;

    sphere.center = {1, 2, 3};
    sphere.radius = 2;

    md::point const pt_inside = {2, 3, 2};
    md::point const pt_outside = {0, 0, 0};
    md::point const pt_surface = {3, 2, 3};

    CHECK(sphere.implicit(pt_inside) < 0);
    CHECK(sphere.implicit(pt_outside) > 0);
    CHECK(sphere.implicit(pt_surface) == Approx(0).margin(1e-6));
}

TEST_CASE("sphere_surface_forcefield - computes inward forcefield")
{
    class inward_forcefield : public md::basic_sphere_surface_forcefield<inward_forcefield>
    {
    public:
        md::harmonic_potential sphere_inward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    md::scalar const x0 = system.add_particle().position.x = -1.2;
    md::scalar const x1 = system.add_particle().position.x = -0.6;
    md::scalar const x2 = system.add_particle().position.x = -0.4;
    md::scalar const x3 = system.add_particle().position.x = 0.4;
    md::scalar const x4 = system.add_particle().position.x = 0.8;
    md::scalar const x5 = system.add_particle().position.x = 1.2;
    md::scalar const x6 = system.add_particle().position.x = 1.5;
    md::scalar const x7 = system.add_particle().position.x = 2.0;

    (void) x0;
    (void) x1;
    (void) x2;
    (void) x7;

    // Particles 2 to 6 are inside this sphere.
    md::sphere sphere;
    sphere.center.x = 0.6;
    sphere.radius = 1.1;

    md::scalar const xl = sphere.center.x - sphere.radius;
    md::scalar const xh = sphere.center.x + sphere.radius;

    inward_forcefield inward;
    inward.set_sphere(sphere);

    // Energy
    md::scalar const e2 = 0.5 * (x2 - xl) * (x2 - xl);
    md::scalar const e3 = 0.5 * (x3 - xl) * (x3 - xl);
    md::scalar const e4 = 0.5 * (x4 - xh) * (x4 - xh);
    md::scalar const e5 = 0.5 * (x5 - xh) * (x5 - xh);
    md::scalar const e6 = 0.5 * (x6 - xh) * (x6 - xh);
    md::scalar const expected_energy = e2 + e3 + e4 + e5 + e6;

    CHECK(inward.compute_energy(system) == Approx(expected_energy));

    // Force
    std::vector<md::vector> forces(system.particle_count());
    inward.compute_force(system, forces);

    CHECK(forces[0].x == 0);
    CHECK(forces[1].x == 0);
    CHECK(forces[2].x == Approx(xl - x2));
    CHECK(forces[3].x == Approx(xl - x3));
    CHECK(forces[4].x == Approx(xh - x4));
    CHECK(forces[5].x == Approx(xh - x5));
    CHECK(forces[6].x == Approx(xh - x6));
    CHECK(forces[7].x == 0);
}

TEST_CASE("sphere_surface_forcefield - computes outward forcefield")
{
    class outward_forcefield : public md::basic_sphere_surface_forcefield<outward_forcefield>
    {
    public:
        md::harmonic_potential sphere_outward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    md::scalar const x0 = system.add_particle().position.x = -1.2;
    md::scalar const x1 = system.add_particle().position.x = -0.6;
    md::scalar const x2 = system.add_particle().position.x = -0.4;
    md::scalar const x3 = system.add_particle().position.x = 0.4;
    md::scalar const x4 = system.add_particle().position.x = 0.8;
    md::scalar const x5 = system.add_particle().position.x = 1.2;
    md::scalar const x6 = system.add_particle().position.x = 1.5;
    md::scalar const x7 = system.add_particle().position.x = 2.0;

    (void) x2;
    (void) x3;
    (void) x4;
    (void) x5;
    (void) x6;

    // Particles 0, 1 and 7 are inside this sphere.
    md::sphere sphere;
    sphere.center.x = 0.6;
    sphere.radius = 1.1;

    md::scalar const xl = sphere.center.x - sphere.radius;
    md::scalar const xh = sphere.center.x + sphere.radius;

    outward_forcefield outward;
    outward.set_sphere(sphere);

    // Energy
    md::scalar const e0 = 0.5 * (x0 - xl) * (x0 - xl);
    md::scalar const e1 = 0.5 * (x1 - xl) * (x1 - xl);
    md::scalar const e7 = 0.5 * (x7 - xh) * (x7 - xh);
    md::scalar const expected_energy = e0 + e1 + e7;

    CHECK(outward.compute_energy(system) == Approx(expected_energy));

    // Force
    std::vector<md::vector> forces(system.particle_count());
    outward.compute_force(system, forces);

    CHECK(forces[0].x == Approx(xl - x0));
    CHECK(forces[1].x == Approx(xl - x1));
    CHECK(forces[2].x == 0);
    CHECK(forces[3].x == 0);
    CHECK(forces[4].x == 0);
    CHECK(forces[5].x == 0);
    CHECK(forces[6].x == 0);
    CHECK(forces[7].x == Approx(xh - x7));
}

TEST_CASE("sphere_surface_forcefield::compute_force - adds force to array")
{
    class outward_forcefield : public md::basic_sphere_surface_forcefield<outward_forcefield>
    {
    public:
        md::harmonic_potential sphere_outward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    // A particle near surface
    system.add_particle().position = {1.1, 0, 0};

    md::sphere sphere;
    sphere.center = {0, 0, 0};
    sphere.radius = 1;

    outward_forcefield outward;
    outward.set_sphere(sphere);

    // compute_force does not clear existing force
    std::vector<md::vector> forces = {
        {1, 2, 3}
    };
    outward.compute_force(system, forces);

    CHECK(forces[0].x == Approx(1 - (1.1 - 1.0)));
    CHECK(forces[0].y == Approx(2));
    CHECK(forces[0].z == Approx(3));
}

TEST_CASE("sphere_surface_forcefield::compute_force - collects normal force stats")
{
    class surface_forcefield : public md::basic_sphere_surface_forcefield<surface_forcefield>
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
        system.add_particle().position = {1.1, 0, 0};
        system.add_particle().position = {0, 1.1, 0};
        system.add_particle().position = {0, 0, 1.1};

        md::sphere sphere;
        sphere.center = {0, 0, 0};
        sphere.radius = 1;

        surface_forcefield ff;
        ff.set_sphere(sphere);

        std::vector<md::vector> forces(system.particle_count());
        ff.compute_force(system, forces);
        auto const sum_normal = forces[0].x + forces[1].y + forces[2].z;

        CHECK(ff.stats.reaction_force == Approx(-sum_normal));
    }

    SECTION("inward force")
    {
        md::system system;
        system.add_particle().position = {0.9, 0, 0};
        system.add_particle().position = {0, 0.9, 0};
        system.add_particle().position = {0, 0, 0.9};

        md::sphere sphere;
        sphere.center = {0, 0, 0};
        sphere.radius = 1;

        surface_forcefield ff;
        ff.set_sphere(sphere);

        std::vector<md::vector> forces(system.particle_count());
        ff.compute_force(system, forces);
        auto const sum_normal = forces[0].x + forces[1].y + forces[2].z;

        CHECK(ff.stats.reaction_force == Approx(-sum_normal));
    }
}

TEST_CASE("make_sphere_inward_forcefield - creates a sphere_surface_forcefield")
{
    SECTION("fixed potential")
    {
        auto ff =
            md::make_sphere_inward_forcefield(
                md::harmonic_potential{1.23}
            )
            .set_sphere(md::sphere{});

        md::system system;
        md::harmonic_potential pot = ff.sphere_inward_potential(system, 0);
        CHECK(pot.spring_constant == 1.23);

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::sphere_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda potential")
    {
        auto ff =
            md::make_sphere_inward_forcefield([](md::index i) {
                return md::harmonic_potential{i * 1.0};
            })
            .set_sphere(md::sphere{});

        md::system system;
        md::harmonic_potential pot1 = ff.sphere_inward_potential(system, 1);
        md::harmonic_potential pot2 = ff.sphere_inward_potential(system, 2);
        CHECK(pot1.spring_constant == Approx(1));
        CHECK(pot2.spring_constant == Approx(2));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::sphere_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda sphere")
    {
        auto ff =
            md::make_sphere_inward_forcefield(
                md::harmonic_potential{1.23}
            )
            .set_sphere([&] {
                return md::sphere{md::point{1, 2, 3}, 4.5};
            });

        md::system system;
        md::sphere sphere = ff.sphere(system);

        CHECK(sphere.center.x == Approx(1));
        CHECK(sphere.center.y == Approx(2));
        CHECK(sphere.center.z == Approx(3));
        CHECK(sphere.radius == Approx(4.5));
    }
}

TEST_CASE("make_sphere_outward_forcefield - creates a sphere_surface_forcefield")
{
    SECTION("fixed potential")
    {
        auto ff =
            md::make_sphere_outward_forcefield(
                md::harmonic_potential{1.23}
            )
            .set_sphere(md::sphere{});

        md::system system;
        md::harmonic_potential pot = ff.sphere_outward_potential(system, 0);
        CHECK(pot.spring_constant == 1.23);

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::sphere_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda potential")
    {
        auto ff =
            md::make_sphere_outward_forcefield([](md::index i) {
                return md::harmonic_potential{i * 1.0};
            })
            .set_sphere(md::sphere{});

        md::system system;
        md::harmonic_potential pot1 = ff.sphere_outward_potential(system, 1);
        md::harmonic_potential pot2 = ff.sphere_outward_potential(system, 2);
        CHECK(pot1.spring_constant == Approx(1));
        CHECK(pot2.spring_constant == Approx(2));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::sphere_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda sphere")
    {
        auto ff =
            md::make_sphere_outward_forcefield(
                md::harmonic_potential{1.23}
            )
            .set_sphere([&] {
                return md::sphere{md::point{1, 2, 3}, 4.5};
            });

        md::system system;
        md::sphere sphere = ff.sphere(system);

        CHECK(sphere.center.x == Approx(1));
        CHECK(sphere.center.y == Approx(2));
        CHECK(sphere.center.z == Approx(3));
        CHECK(sphere.radius == Approx(4.5));
    }
}
