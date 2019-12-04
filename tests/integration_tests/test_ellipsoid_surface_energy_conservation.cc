#include <md/basic_types.hpp>
#include <md/system.hpp>
#include <md/forcefield/ellipsoid_surface_forcefield.hpp>
#include <md/potential/harmonic_potential.hpp>
#include <md/simulation/newtonian_dynamics.hpp>

#include <catch.hpp>


TEST_CASE("Energy conservation of particles in a ellipsoid")
{
    class ellipsoid_packing_forcefield
        : public md::basic_ellipsoid_surface_forcefield<ellipsoid_packing_forcefield>
    {
    public:
        auto ellipsoid_outward_potential(md::system const&, md::index)
        {
            md::harmonic_potential harmonic;
            harmonic.spring_constant = 1000;
            return harmonic;
        }
    };

    md::system system;

    md::ellipsoid ellipsoid;
    ellipsoid.semiaxis_x = 0.6;
    ellipsoid.semiaxis_y = 0.4;
    ellipsoid.semiaxis_z = 0.2;

    md::particle_ref part = system.add_particle();
    part.mass = 0.1;
    part.position = ellipsoid.center;
    part.velocity = {1, 2, 3};

    ellipsoid_packing_forcefield forcefield;
    forcefield.set_ellipsoid(ellipsoid);
    system.add_forcefield(forcefield);

    md::newtonian_dynamics_config config;
    config.timestep = 1e-4;
    config.steps = 100000;

    md::scalar const initial_energy = system.compute_energy();
    md::simulate_newtonian_dynamics(system, config);
    md::scalar const final_energy = system.compute_energy();

    CHECK(final_energy == Approx(initial_energy).epsilon(0.001));
}
