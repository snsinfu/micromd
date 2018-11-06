#include <functional>
#include <memory>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/simulation/brownian_dynamics.hpp>

#include <catch.hpp>


TEST_CASE("brownian_dynamics_config - defaults to sane configuration")
{
    md::brownian_dynamics_config config;

    CHECK(config.temperature == 1);
    CHECK(config.timestep == 1);
    CHECK(config.steps == 1);
    CHECK(config.seed == 0);
    CHECK(!config.callback);
}

TEST_CASE("simulate_brownian_dynamics - does nothing if steps is zero")
{
    md::system system;

    system.add_particle();

    md::brownian_dynamics_config config;
    config.steps = 0;

    md::simulate_brownian_dynamics(system, config);

    CHECK(system.view_positions()[0].x == 0);
    CHECK(system.view_positions()[0].y == 0);
    CHECK(system.view_positions()[0].z == 0);
}

TEST_CASE("simulate_brownian_dynamics - can simulate frozen system")
{
    md::system system;

    system.add_particle();

    md::brownian_dynamics_config config;
    config.temperature = 0;
    config.steps = 100;

    md::simulate_brownian_dynamics(system, config);

    CHECK(system.view_positions()[0].x == 0);
    CHECK(system.view_positions()[0].y == 0);
    CHECK(system.view_positions()[0].z == 0);
}

TEST_CASE("simulate_brownian_dynamics - can simulate brownian motion")
{
    md::system system;

    system.add_particle();

    md::brownian_dynamics_config config;
    config.temperature = 1;
    config.steps = 100;

    md::simulate_brownian_dynamics(system, config);

    CHECK(system.view_positions()[0].x != 0);
    CHECK(system.view_positions()[0].y != 0);
    CHECK(system.view_positions()[0].z != 0);
}

TEST_CASE("simulate_brownian_dynamics - supports adaptive time-stepping")
{
    struct motion_recorder
    {
        md::system* system;
        md::scalar sum_motion = 0;
        md::step steps = 0;
        md::point prev_position;

        void operator()(md::step)
        {
            md::point const position = system->view_positions()[0];
            sum_motion += (position - prev_position).norm();
            prev_position = position;
            steps++;
        }

        md::scalar mean() const
        {
            return sum_motion / md::scalar(steps);
        }
    };

    md::system system;
    system.add_particle();
    system.add_particle();
    system.add_particle();

    motion_recorder recorder;
    recorder.system = &system;

    md::brownian_dynamics_config config;
    config.temperature = 1;
    config.spacestep = 0.01;
    config.steps = 1000;
    config.callback = std::ref(recorder);

    md::simulate_brownian_dynamics(system, config);

    CHECK(recorder.mean() == Approx(config.spacestep).epsilon(0.1));
}
