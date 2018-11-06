#include <algorithm>
#include <cassert>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/forcefield/sphere_surface_forcefield.hpp>
#include <md/potential/harmonic_potential.hpp>

#include <catch.hpp>


TEST_CASE("sphere_surface_forcefield - can force inward surface")
{
    class inward_forcefield : public md::sphere_surface_forcefield<inward_forcefield>
    {
    public:
        md::sphere sphere_surface(md::system const&)
        {
            return {};
        }

        md::harmonic_potential sphere_inward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    inward_forcefield inward;

    md::system system;

    system.add_particle().position = {-1.2, 0, 0};
    system.add_particle().position = {-0.4, 0, 0};
    system.add_particle().position = { 0.4, 0, 0};
    system.add_particle().position = { 1.2, 0, 0};

    inward.compute_energy(system);
}
