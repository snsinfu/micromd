#include <md/basic_types.hpp>
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
