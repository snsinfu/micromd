#include <algorithm>
#include <vector>

#include <md/forcefield.hpp>
#include <md/system.hpp>
#include <md/typedef.hpp>

#include <md/forcefield/neighbor_pair_forcefield.hpp>

#include "../catch.hpp"


namespace
{
    struct test_potential
    {
        md::scalar evaluate_energy(md::vector r) const
        {
            if (r.squared_norm() < 1) {
                return 1 - r.squared_norm();
            } else {
                return 0;
            }
        }

        md::vector evaluate_force(md::vector r) const
        {
            if (r.squared_norm() < 1) {
                return 2 * r;
            } else {
                return {};
            }
        }
    };
}


TEST_CASE("neighbor_pair_forcefield - computes only neighbor interactions")
{
    class my_forcefield : public md::neighbor_pair_forcefield<my_forcefield>
    {
    public:
        md::scalar neighbor_distance(md::system const&)
        {
            return 1;
        }

        test_potential neighbor_pair_potential(md::system const&, md::index, md::index)
        {
            return test_potential{};
        }
    };

    // Points on a line. Only adjacent points are relevant.
    md::system system;

    system.add_particle();
    system.add_particle();
    system.add_particle();
    system.add_particle();

    md::array_view<md::point> positions = system.view(md::position_attribute);

    positions[0] = {0.0, 0, 0};
    positions[1] = {0.6, 0, 0};
    positions[2] = {1.2, 0, 0};
    positions[3] = {1.8, 0, 0};

    SECTION("energy is zero")
    {
        my_forcefield forcefield;
        test_potential potential;

        md::scalar const expected_interaction = potential.evaluate_energy({0.6, 0, 0});

        // There are 3 interactions
        CHECK(forcefield.compute_energy(system) == Approx(3 * expected_interaction));
    }

    SECTION("force is zero")
    {
        my_forcefield forcefield;

        std::vector<md::vector> forces(system.particle_count());
        forcefield.compute_force(system, forces);

        CHECK(forces[0].x == Approx(-1.2));
        CHECK(forces[1].x == Approx(0).margin(1e-6));
        CHECK(forces[2].x == Approx(0).margin(1e-6));
        CHECK(forces[3].x == Approx(1.2));
    }
}
