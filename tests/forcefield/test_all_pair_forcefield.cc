#include <algorithm>
#include <cassert>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/forcefield/all_pair_forcefield.hpp>

#include <catch.hpp>


namespace
{
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

TEST_CASE("all_pair_forcefield - computes correct forcefield")
{
    // u(r) = r^2
    struct harmonic_potential
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

    class test_forcefield : public md::all_pair_forcefield<test_forcefield>
    {
    public:
        harmonic_potential all_pair_potential(md::system const&, md::index, md::index)
        {
            return harmonic_potential{};
        }
    };

    // Particles on a 5x5x5 grid
    md::system system;

    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            for (int z = -2; z <= 2; z++) {
                system.add_particle().position = {
                    x * 0.5,
                    y * 0.5,
                    z * 0.5
                };
            }
        }
    }

    test_forcefield forcefield;
    harmonic_potential potential;

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
