#include <algorithm>
#include <cassert>
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

    md::scalar max_difference(
        md::array_view<md::vector const> v1,
        md::array_view<md::vector const> v2
    )
    {
        assert(v1.size() == v2.size());

        md::scalar diff = 0;

        for (md::index i = 0; i < v1.size(); i++) {
            diff = std::max(diff, md::norm(v1[i] - v2[i]));
        }
        return diff;
    }
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

    test_forcefield forcefield;
    bell_potential potential;

    std::vector<md::vector> actual_forces(system.particle_count());
    std::vector<md::vector> expected_forces(system.particle_count());

    forcefield.compute_force(system, actual_forces);

    md::scalar actual_energy = forcefield.compute_energy(system);
    md::scalar expected_energy = 0;

    md::array_view<md::point const> positions = system.view_positions();

    for (md::index i = 0; i < positions.size(); i++) {
        for (md::index j = 0; j < i; j++) {
            md::vector const r = positions[i] - positions[j];
            md::vector const force = potential.evaluate_force(r);
            md::scalar const energy = potential.evaluate_energy(r);

            expected_forces[i] += force;
            expected_forces[j] -= force;

            expected_energy += energy;
        }
    }

    CHECK(actual_energy == Approx(expected_energy));
    CHECK(max_difference(actual_forces, expected_forces) < 1e-6);
}

TEST_CASE("neighbor_pair_forcefield::compute_force - adds force to array")
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

    md::system system;

    system.add_particle().position = {0.5, 0, 0};
    system.add_particle().position = {0, 0.5, 0};

    test_forcefield ff;

    // compute_force does not clear existing force
    std::vector<md::vector> forces = {
        {1, 2, 3},
        {4, 5, 6}
    };
    ff.compute_force(system, forces);

    CHECK(forces[0].x == Approx(1 + 2 * 0.5));
    CHECK(forces[0].y == Approx(2 - 2 * 0.5));
    CHECK(forces[0].z == Approx(3));
    CHECK(forces[1].x == Approx(4 - 2 * 0.5));
    CHECK(forces[1].y == Approx(5 + 2 * 0.5));
    CHECK(forces[1].z == Approx(6));
}

TEST_CASE("make_neighbor_pair_forcefield - creates a neighbor_pair_forcefield")
{
    struct harmonic_potential
    {
        md::scalar spring_constant;

        md::scalar evaluate_energy(md::vector r) const
        {
            return spring_constant * r.squared_norm() / 2;
        }

        md::vector evaluate_force(md::vector r) const
        {
            return -spring_constant * r;
        }
    };

    md::system system;

    auto ff = md::make_neighbor_pair_forcefield(harmonic_potential{1.23});
    ff.set_neighbor_distance(4.56);

    auto ndist = ff.neighbor_distance(system);
    auto pot = ff.neighbor_pair_potential(system, 0, 1);

    using ff_type = decltype(ff);
    using pot_type = decltype(pot);

    CHECK(std::is_base_of<md::neighbor_pair_forcefield<ff_type>, ff_type>::value);
    CHECK(std::is_same<pot_type, harmonic_potential>::value);
    CHECK(pot.spring_constant == 1.23);
    CHECK(ndist == 4.56);
}
