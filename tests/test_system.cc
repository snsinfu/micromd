#include "../include/md.hpp"
#include "catch.hpp"


TEST_CASE("system - is default constructible")
{
    md::system system;

    CHECK(system.particle_count() == 0);
}

TEST_CASE("system::add_particle - adds a particle")
{
    md::system system;

    system.add_particle();
    CHECK(system.particle_count() == 1);

    system.add_particle();
    CHECK(system.particle_count() == 2);
}
