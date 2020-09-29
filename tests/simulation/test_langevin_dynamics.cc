#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/simulation/langevin_dynamics.hpp>

#include <catch.hpp>


TEST_CASE("langevin_dynamics_config - defaults to sane configuration")
{
    md::langevin_dynamics_config config;

    CHECK(config.temperature == 1);
    CHECK(config.timestep == 1);
    CHECK(config.steps == 1);
    CHECK(config.seed == 0);
    CHECK(!config.callback);
}

TEST_CASE("simulate_langevin_dynamics - callback step is 1-based")
{
    md::system system;
    system.add_particle();

    std::vector<md::step> actual_steps;
    std::vector<md::step> expected_steps = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    md::langevin_dynamics_config config;
    config.steps = 10;
    config.callback = [&](md::step step) { actual_steps.push_back(step); };

    md::simulate_langevin_dynamics(system, config);

    CHECK(actual_steps == expected_steps);
}

TEST_CASE("simulate_langevin_dynamics - modifies particle velocity")
{
    md::system system;
    system.add_particle();

    md::langevin_dynamics_config config;
    config.steps = 1000;
    md::simulate_langevin_dynamics(system, config);

    md::vector const velocity = system.view_velocities()[0];
    CHECK(velocity.norm() > 0);
}

TEST_CASE("simulate_langevin_dynamics - produces equipartition")
{
    // Simulate athermal dynamics of single particle in a harmonic well. We
    // sample kinetic/potential energy and test equipartition property.

    md::system system;
    {
        md::particle_ref part = system.add_particle();
        part.mass = 1.23;
        part.friction = 3.21e4;
    }

    class harmonic_forcefield : public md::forcefield
    {
    public:
        explicit
        harmonic_forcefield(md::scalar spring_constant)
            : spring_constant_{spring_constant}
        {
        }

        md::scalar
        compute_energy(md::system const& system) override
        {
            md::array_view<md::point const> positions = system.view_positions();
            md::point const pos = positions[0];
            return 0.5 * spring_constant_ * pos.vector().squared_norm();
        }

        void
        compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();
            md::point const pos = positions[0];
            md::vector const force = -spring_constant_ * pos.vector();
            forces[0] += force;
        }

    private:
        md::scalar spring_constant_;
    };

    system.add_forcefield(harmonic_forcefield(7.89));

    md::langevin_dynamics_config config;
    config.temperature = 1.11;
    config.timestep = 0.01;
    config.steps = 100000;
    config.seed = 2;

    md::scalar mean_kinetic_energy = 0;
    md::scalar mean_potential_energy = 0;
    config.callback = [&](md::step) {
        mean_kinetic_energy += system.compute_kinetic_energy();
        mean_potential_energy += system.compute_potential_energy();
    };

    md::simulate_langevin_dynamics(system, config);

    mean_kinetic_energy /= md::scalar(config.steps);
    mean_potential_energy /= md::scalar(config.steps);

    // The equipartition theorem
    CHECK(mean_kinetic_energy == Approx(1.5 * config.temperature).epsilon(0.1));
    CHECK(mean_potential_energy == Approx(1.5 * config.temperature).epsilon(0.1));
}
