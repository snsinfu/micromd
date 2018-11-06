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
    class surface_forcefield : public md::sphere_surface_forcefield<surface_forcefield>
    {
    public:
        md::harmonic_potential sphere_inward_potential(md::system const&, md::index)
        {
            return md::harmonic_potential{};
        }
    };

    md::sphere sphere;
    sphere.center = {0, 0, 1.1};
    sphere.radius = 1.1;

    surface_forcefield surface;
    surface.set_sphere(sphere);

    md::system system;

    system.add_particle().position = {0, 0, -1.2};
    system.add_particle().position = {0, 0, -0.4};
    system.add_particle().position = {0, 0,  0.4};
    system.add_particle().position = {0, 0,  1.2};

    surface.compute_energy(system);
}
