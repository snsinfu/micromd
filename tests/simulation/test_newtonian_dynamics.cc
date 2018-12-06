#include <memory>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/simulation/newtonian_dynamics.hpp>

#include <catch.hpp>


namespace
{
    class gravitational_forcefield : public md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const& system) override
        {
            md::scalar sum = 0;
            for (md::point pos : system.view_positions()) {
                sum += pos.z;
            }
            return sum;
        }

        void compute_force(md::system const&, md::array_view<md::vector> forces) override
        {
            for (md::vector& force : forces) {
                force.z -= 1;
            }
        }
    };
}


TEST_CASE("newtonian_dynamics_config - defaults to sane configuration")
{
    md::newtonian_dynamics_config config;

    CHECK(config.timestep == 1);
    CHECK(config.steps == 1);
    CHECK(!config.callback);
}

TEST_CASE("simulate_newtonian_dynamics - does nothing if steps is zero")
{
    md::system system;
    system.add_particle();
    system.add_forcefield(std::make_shared<::gravitational_forcefield>());

    md::newtonian_dynamics_config config;
    config.steps = 0;

    md::simulate_newtonian_dynamics(system, config);

    CHECK(system.view_positions()[0].x == 0);
    CHECK(system.view_positions()[0].y == 0);
    CHECK(system.view_positions()[0].z == 0);
}

TEST_CASE("simulate_newtonian_dynamics - calls callback function on each step")
{
    md::system system;
    system.add_particle();

    int counter = 0;

    md::newtonian_dynamics_config config;
    config.steps = 10;
    config.callback = [&](md::step) { counter++; };

    md::simulate_newtonian_dynamics(system, config);

    CHECK(counter == 10);
}

TEST_CASE("simulate_newtonian_dynamics - callback step is 1-based")
{
    md::system system;
    system.add_particle();

    std::vector<md::step> actual_steps;
    std::vector<md::step> expected_steps = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    md::newtonian_dynamics_config config;
    config.steps = 10;
    config.callback = [&](md::step step) { actual_steps.push_back(step); };

    md::simulate_newtonian_dynamics(system, config);

    CHECK(actual_steps == expected_steps);
}
