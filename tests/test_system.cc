#include "../include/md.hpp"
#include "catch.hpp"


TEST_CASE("system::particle_count - returns zero on default state")
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

TEST_CASE("system::add_particle - extends existing attribute array")
{
    md::system system;

    // Add attributes
    md::attribute_key<int, struct test_attribute_a_tag> test_attribute_a = {};
    md::attribute_key<int, struct test_attribute_b_tag> test_attribute_b = {};

    system.require(test_attribute_a);
    system.require(test_attribute_b);

    CHECK(system.view(test_attribute_a).size() == 0);
    CHECK(system.view(test_attribute_b).size() == 0);

    // Extend by one
    system.add_particle();

    CHECK(system.view(test_attribute_a).size() == 1);
    CHECK(system.view(test_attribute_b).size() == 1);

    // More
    system.add_particle();
    system.add_particle();

    CHECK(system.view(test_attribute_a).size() == 3);
    CHECK(system.view(test_attribute_b).size() == 3);
}

TEST_CASE("system::require - creates an attribute if it does not exist")
{
    // Add 3 particles
    md::system system;

    system.add_particle();
    system.add_particle();
    system.add_particle();

    md::attribute_key<int, struct test_attribute_tag> test_attribute = {};

    // This call creates an attribute with 3 entries
    system.require(test_attribute);

    md::array_view<int> values = system.view(test_attribute);
    CHECK(values.size() == system.particle_count());
    CHECK(values[0] == 0);
    CHECK(values[1] == 0);
    CHECK(values[2] == 0);

    // Update attribute values and check that an extra require() call does not
    // tamper with the existing attribute.
    values[0] = 42;
    values[1] = 42;
    values[2] = 42;

    system.require(test_attribute);

    CHECK(values[0] == 42);
    CHECK(values[1] == 42);
    CHECK(values[2] == 42);
}
