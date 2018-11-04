#include <algorithm>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/forcefield/neighbor_pair_forcefield.hpp>

#include <catch.hpp>


namespace
{
    // u(r) = 1 - r^2  (r < 1)
    struct bell_potential
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


TEST_CASE("neighbor_pair_forcefield - computes correct forcefield")
{
    class test_forcefield : public md::neighbor_pair_forcefield<test_forcefield>
    {
    public:
        md::scalar neighbor_distance(md::system const&)
        {
            return 1;
        }

        bell_potential neighbor_pair_potential(md::system const&, md::index, md::index)
        {
            return bell_potential{};
        }
    };

    // Particles on a 5x5x5 grid
    md::system system;

    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            for (int z = -2; z <= 2; z++) {
                md::basic_particle_data data;
                data.position = {
                    x * 0.5,
                    y * 0.5,
                    z * 0.5
                };
                system.add_particle(data);
            }
        }
    }

    SECTION("energy is correct")
    {
        bell_potential potential;
        test_forcefield forcefield;

        md::scalar expected = 0;
        md::array_view<md::point const> positions = system.view_positions();

        for (md::index i = 0; i < positions.size(); i++) {
            for (md::index j = 0; j < i; j++) {
                expected += potential.evaluate_energy(positions[i] - positions[j]);
            }
        }

        CHECK(forcefield.compute_energy(system) == expected);
    }

    SECTION("force is correct")
    {
        // FIXME
    }
}

TEST_CASE("force_neighbor_pairs")
{
    // Particles on a 5x5x5 grid
    md::system system;

    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            for (int z = -2; z <= 2; z++) {
                md::basic_particle_data data;
                data.position = {
                    x * 0.5,
                    y * 0.5,
                    z * 0.5
                };
                system.add_particle(data);
            }
        }
    }

    SECTION("energy is correct")
    {
        md::force_neighbor_pairs(system, 1, [](md::index, md::index) {
            return bell_potential{};
        });
        bell_potential potential;

        md::scalar expected = 0;
        md::array_view<md::point const> positions = system.view_positions();

        for (md::index i = 0; i < positions.size(); i++) {
            for (md::index j = 0; j < i; j++) {
                expected += potential.evaluate_energy(positions[i] - positions[j]);
            }
        }

        CHECK(system.compute_potential_energy() == expected);
    }

    SECTION("force is correct")
    {
        // FIXME
    }
}
