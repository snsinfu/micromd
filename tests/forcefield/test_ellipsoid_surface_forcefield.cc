#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <md/forcefield/ellipsoid_surface_forcefield.hpp>

#include <catch.hpp>


TEST_CASE("ellipsoid::implicit - returns negative inside, positive outside")
{
    md::ellipsoid ellip;

    ellip.center = {1, 2, 3};
    ellip.semiaxis_x = 4;
    ellip.semiaxis_y = 5;
    ellip.semiaxis_z = 6;

    md::point const pt_inside = {3, 3, 3};
    md::point const pt_outside = {-3, -3, -3};
    md::point const pt_surface = {1, 7, 3};

    CHECK(ellip.implicit(pt_inside) < 0);
    CHECK(ellip.implicit(pt_outside) > 0);
    CHECK(ellip.implicit(pt_surface) == Approx(0).margin(1e-6));
}

TEST_CASE("detail::evaluate_point - evaluates a point in an ellipsoid")
{
    md::ellipsoid ellip;

    SECTION("center is undefined")
    {
        md::detail::ellipsoid_eval const ev = md::detail::evaluate_point(ellip, ellip.center);

        CHECK(ev.undefined);
    }

    SECTION("surface")
    {
        md::point const pt = {1, 0, 0};
        md::detail::ellipsoid_eval const ev = md::detail::evaluate_point(ellip, pt);

        CHECK_FALSE(ev.undefined);

        CHECK(ev.delta.x == 0);
        CHECK(ev.delta.y == 0);
        CHECK(ev.delta.z == 0);

        CHECK(ev.strain.x == 0);
        CHECK(ev.strain.y == 0);
        CHECK(ev.strain.z == 0);

        CHECK(ev.implicit == 0);
    }

    SECTION("inner skin")
    {
        md::point const pt = {0.99, 0, 0};
        md::detail::ellipsoid_eval const ev = md::detail::evaluate_point(ellip, pt);

        CHECK_FALSE(ev.undefined);

        CHECK(ev.delta.x == Approx(-0.01).epsilon(0.1));
        CHECK(ev.delta.y == 0);
        CHECK(ev.delta.z == 0);

        CHECK(ev.strain.x == Approx(-0.0101520));
        CHECK(ev.strain.y == Approx(-0.0101520));
        CHECK(ev.strain.z == Approx(-0.0101520));

        CHECK(ev.implicit == Approx(-0.0199));
    }

    SECTION("outer skin")
    {
        md::point const pt = {1.01, 0, 0};
        md::detail::ellipsoid_eval const ev = md::detail::evaluate_point(ellip, pt);

        CHECK_FALSE(ev.undefined);

        CHECK(ev.delta.x == Approx(0.01).epsilon(0.1));
        CHECK(ev.delta.y == 0);
        CHECK(ev.delta.z == 0);

        CHECK(ev.strain.x == Approx(0.00985198));
        CHECK(ev.strain.y == Approx(0.00985198));
        CHECK(ev.strain.z == Approx(0.00985198));

        CHECK(ev.implicit == Approx(0.0201));
    }
}

TEST_CASE("ellipsoid_surface_forcefield - computes inward forcefield")
{
    class inward_forcefield : public md::ellipsoid_surface_forcefield<inward_forcefield>
    {
    public:
        md::harmonic_potential ellipsoid_inward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::ellipsoid ellipsoid;
    ellipsoid.center = {0.5, 0, 0};
    ellipsoid.semiaxis_x = 1;
    ellipsoid.semiaxis_y = 2;
    ellipsoid.semiaxis_z = 3;

    md::scalar const xl = ellipsoid.center.x - ellipsoid.semiaxis_x;
    md::scalar const xh = ellipsoid.center.x + ellipsoid.semiaxis_x;

    md::system system;
    md::scalar const x0 = system.add_particle().position.x = -0.60;
    md::scalar const x1 = system.add_particle().position.x = -0.49;
    md::scalar const x2 = system.add_particle().position.x = 0.50;
    md::scalar const x3 = system.add_particle().position.x = 1.49;
    md::scalar const x4 = system.add_particle().position.x = 1.60;

    (void) x0;
    (void) x2; // No interaction is calculated for x3 == center
    (void) x4;

    inward_forcefield inward;
    inward.set_ellipsoid(ellipsoid);

    // Energy
    md::scalar const e1 = 0.5 * (xl - x1) * (xl - x1);
    md::scalar const e3 = 0.5 * (xh - x3) * (xh - x3);
    md::scalar const expected_energy = e1 + e3;

    CHECK(inward.compute_energy(system) == Approx(expected_energy).epsilon(0.1));

    // Force
    std::vector<md::vector> forces(system.particle_count());

    md::scalar const f1 = xl - x1;
    md::scalar const f3 = xh - x3;

    inward.compute_force(system, forces);

    CHECK(forces[0].x == 0);
    CHECK(forces[1].x == Approx(f1).epsilon(0.1));
    CHECK(forces[2].x == 0);
    CHECK(forces[3].x == Approx(f3).epsilon(0.1));
    CHECK(forces[4].x == 0);
}

TEST_CASE("ellipsoid_surface_forcefield - computes outward forcefield")
{
    class outward_forcefield : public md::ellipsoid_surface_forcefield<outward_forcefield>
    {
    public:
        md::harmonic_potential ellipsoid_outward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::ellipsoid ellipsoid;
    ellipsoid.center = {0.5, 0, 0};
    ellipsoid.semiaxis_x = 1;
    ellipsoid.semiaxis_y = 2;
    ellipsoid.semiaxis_z = 3;

    md::scalar const xl = ellipsoid.center.x - ellipsoid.semiaxis_x;
    md::scalar const xh = ellipsoid.center.x + ellipsoid.semiaxis_x;

    md::system system;
    md::scalar const x0 = system.add_particle().position.x = -0.51;
    md::scalar const x1 = system.add_particle().position.x = -0.49;
    md::scalar const x2 = system.add_particle().position.x = 0.50;
    md::scalar const x3 = system.add_particle().position.x = 1.49;
    md::scalar const x4 = system.add_particle().position.x = 1.51;

    (void) x1;
    (void) x2;
    (void) x3;

    outward_forcefield outward;
    outward.set_ellipsoid(ellipsoid);

    // Energy
    md::scalar const e0 = 0.5 * (xl - x0) * (xl - x0);
    md::scalar const e4 = 0.5 * (xh - x4) * (xh - x4);
    md::scalar const expected_energy = e0 + e4;

    CHECK(outward.compute_energy(system) == Approx(expected_energy).epsilon(0.1));

    // Force
    std::vector<md::vector> forces(system.particle_count());

    md::scalar const f0 = xl - x0;
    md::scalar const f4 = xh - x4;

    outward.compute_force(system, forces);

    CHECK(forces[0].x == Approx(f0).epsilon(0.1));
    CHECK(forces[1].x == 0);
    CHECK(forces[2].x == 0);
    CHECK(forces[3].x == 0);
    CHECK(forces[4].x == Approx(f4).epsilon(0.1));
}

TEST_CASE("ellipsoid_surface_forcefield::set_ellipsoid - returns self")
{
    class test_forcefield : public md::ellipsoid_surface_forcefield<test_forcefield>
    {
    public:
        md::harmonic_potential ellipsoid_outward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    test_forcefield test;
    test_forcefield& ref = test.set_ellipsoid(md::ellipsoid{});

    CHECK(&ref == &test);
}

TEST_CASE("ellipsoid_surface_forcefield::compute_force - adds force to array")
{
    class outward_forcefield : public md::ellipsoid_surface_forcefield<outward_forcefield>
    {
    public:
        md::harmonic_potential ellipsoid_outward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::system system;

    // A particle near surface
    system.add_particle().position = {1.01, 0, 0};

    md::ellipsoid ellip;
    ellip.center = {0, 0, 0};
    ellip.semiaxis_x = 1;
    ellip.semiaxis_y = 0.8;
    ellip.semiaxis_z = 0.6;

    outward_forcefield outward;
    outward.set_ellipsoid(ellip);

    // compute_force does not clear existing force
    std::vector<md::vector> forces = {
        {1, 2, 3}
    };
    outward.compute_force(system, forces);

    CHECK(forces[0].x == Approx(1 - (1.01 - 1.0)).epsilon(0.1));
    CHECK(forces[0].y == Approx(2).epsilon(0.1));
    CHECK(forces[0].z == Approx(3).epsilon(0.1));
}

TEST_CASE("ellipsoid_surface_forcefield::compute_force - collects normal force stats")
{
    class surface_forcefield : public md::ellipsoid_surface_forcefield<surface_forcefield>
    {
    public:
        md::harmonic_potential ellipsoid_inward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }

        md::harmonic_potential ellipsoid_outward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    SECTION("outward force")
    {
        md::system system;
        system.add_particle().position = {1.1, 0, 0};
        system.add_particle().position = {0, 0.9, 0};
        system.add_particle().position = {0, 0, 0.7};

        md::ellipsoid ellip;
        ellip.center = {0, 0, 0};
        ellip.semiaxis_x = 1;
        ellip.semiaxis_y = 0.8;
        ellip.semiaxis_z = 0.6;

        surface_forcefield ff;
        ff.set_ellipsoid(ellip);

        std::vector<md::vector> forces(system.particle_count());
        ff.compute_force(system, forces);
        auto const sum_normal = forces[0].x + forces[1].y + forces[2].z;

        CHECK(ff.stats.reaction_force == Approx(-sum_normal));
    }

    SECTION("inward force")
    {
        md::system system;
        system.add_particle().position = {0.9, 0, 0};
        system.add_particle().position = {0, 0.7, 0};
        system.add_particle().position = {0, 0, 0.5};

        md::ellipsoid ellip;
        ellip.center = {0, 0, 0};
        ellip.semiaxis_x = 1;
        ellip.semiaxis_y = 0.8;
        ellip.semiaxis_z = 0.6;

        surface_forcefield ff;
        ff.set_ellipsoid(ellip);

        std::vector<md::vector> forces(system.particle_count());
        ff.compute_force(system, forces);
        auto const sum_normal = forces[0].x + forces[1].y + forces[2].z;

        CHECK(ff.stats.reaction_force == Approx(-sum_normal));
    }
}

TEST_CASE("make_ellipsoid_inward_forcefield - creates a ellipsoid_surface_forcefield")
{
    SECTION("fixed potential")
    {
        auto ff = md::make_ellipsoid_inward_forcefield(md::harmonic_potential{1.23});
        ff.set_ellipsoid(md::ellipsoid{});

        md::system system;
        md::harmonic_potential pot = ff.ellipsoid_inward_potential(system, 0);
        CHECK(pot.spring_constant == 1.23);

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::ellipsoid_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda potential")
    {
        auto ff = md::make_ellipsoid_inward_forcefield([](md::index i) {
            return md::harmonic_potential{i * 1.0};
        });
        ff.set_ellipsoid(md::ellipsoid{});

        md::system system;
        md::harmonic_potential pot1 = ff.ellipsoid_inward_potential(system, 1);
        md::harmonic_potential pot2 = ff.ellipsoid_inward_potential(system, 2);
        CHECK(pot1.spring_constant == Approx(1));
        CHECK(pot2.spring_constant == Approx(2));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::ellipsoid_surface_forcefield<ff_type>, ff_type>::value);
    }
}

TEST_CASE("make_ellipsoid_outward_forcefield - creates a ellipsoid_surface_forcefield")
{
    SECTION("fixed potential")
    {
        auto ff = md::make_ellipsoid_outward_forcefield(md::harmonic_potential{1.23});
        ff.set_ellipsoid(md::ellipsoid{});

        md::system system;
        md::harmonic_potential pot = ff.ellipsoid_outward_potential(system, 0);
        CHECK(pot.spring_constant == 1.23);

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::ellipsoid_surface_forcefield<ff_type>, ff_type>::value);
    }

    SECTION("lambda potentialo")
    {
        auto ff = md::make_ellipsoid_outward_forcefield([](md::index i) {
            return md::harmonic_potential{i * 1.0};
        });
        ff.set_ellipsoid(md::ellipsoid{});

        md::system system;
        md::harmonic_potential pot1 = ff.ellipsoid_outward_potential(system, 1);
        md::harmonic_potential pot2 = ff.ellipsoid_outward_potential(system, 2);
        CHECK(pot1.spring_constant == Approx(1));
        CHECK(pot2.spring_constant == Approx(2));

        using ff_type = decltype(ff);
        CHECK(std::is_base_of<md::ellipsoid_surface_forcefield<ff_type>, ff_type>::value);
    }
}
