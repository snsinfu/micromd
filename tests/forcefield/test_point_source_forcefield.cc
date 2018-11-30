#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <md/forcefield/point_source_forcefield.hpp>

#include <catch.hpp>


TEST_CASE("point_source_forcefield - computes correct forcefield")
{
    struct test_potential
    {
        md::scalar evaluate_energy(md::vector r) const
        {
            return r.squared_norm();
        }

        md::vector evaluate_force(md::vector r) const
        {
            return -r;
        }
    };

    class test_forcefield : public md::point_source_forcefield<test_forcefield>
    {
    public:
        test_potential point_source_potential(md::system const&, md::index)
        {
            return test_potential{};
        }
    };

    md::system system;
    test_forcefield forcefield;

    md::point const p0 = system.add_particle().position = {0, 1, 2};
    md::point const p1 = system.add_particle().position = {3, 4, 5};
    md::point const p2 = system.add_particle().position = {6, 7, 8};
    md::point const ps = {9, 0, 1};
    forcefield.set_point_source(ps);

    md::scalar const expected_energy =
        (p0 - ps).squared_norm() +
        (p1 - ps).squared_norm() +
        (p2 - ps).squared_norm();

    CHECK(forcefield.compute_energy(system) == Approx(expected_energy));

    // Force
    std::vector<md::vector> forces(system.particle_count());
    forcefield.compute_force(system, forces);

    CHECK((forces[0] - (ps - p0)).norm() == Approx(0).margin(1e-6));
    CHECK((forces[1] - (ps - p1)).norm() == Approx(0).margin(1e-6));
    CHECK((forces[2] - (ps - p2)).norm() == Approx(0).margin(1e-6));
}
