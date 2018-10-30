#include <memory>

#include "../include/md.hpp"
#include "catch.hpp"


TEST_CASE("system::time - returns zero on default state")
{
    md::system system;

    CHECK(system.time() == 0);
}

TEST_CASE("system::advance_time - adds delta to system time")
{
    md::system system;

    system.advance_time(1.23);
    CHECK(system.time() == Approx(1.23));

    system.advance_time(4.56);
    CHECK(system.time() == Approx(1.23 + 4.56));
}

TEST_CASE("system::advance_time - accepts negative delta")
{
    md::system system;

    system.advance_time(-1.23);
    CHECK(system.time() == Approx(-1.23));

    system.advance_time(4.56);
    CHECK(system.time() == Approx(-1.23 + 4.56));

    system.advance_time(-7.89);
    CHECK(system.time() == Approx(-1.23 + 4.56 - 7.89));
}

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

TEST_CASE("system::add_forcefield - accepts a shared_ptr of a forcefield")
{
    class my_forcefield : public md::forcefield
    {
        md::scalar compute_energy(md::system const&) override
        {
            return 0;
        }

        void compute_force(md::system const&, md::array_view<md::vector>) override
        {
        }
    };

    md::system system;
    system.add_forcefield(std::make_shared<my_forcefield>());
}

TEST_CASE("system::compute_potential_energy - returns zero if no force field is added")
{
    md::system system;

    CHECK(system.compute_potential_energy() == 0);
}

TEST_CASE("system::compute_force - returns zero if no force field is added")
{
    md::system system;
    system.add_particle();

    std::vector<md::vector> forces(system.particle_count(), md::vector{1, 1, 1});
    system.compute_force(forces);

    CHECK(forces[0].x == 0);
    CHECK(forces[0].y == 0);
    CHECK(forces[0].z == 0);
}

TEST_CASE("system::compute_potential_energy - returns the sum of component energy values")
{
    class my_forcefield : public md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const&) override
        {
            return 1;
        }

        void compute_force(md::system const&, md::array_view<md::vector>) override
        {
        }
    };

    md::system system;

    system.add_forcefield(std::make_shared<my_forcefield>());
    system.add_forcefield(std::make_shared<my_forcefield>());

    CHECK(system.compute_potential_energy() == 2);
}

TEST_CASE("system::compute_force - returns the sum of component force vectors")
{
    class my_forcefield : public md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const&) override
        {
            return 1;
        }

        void compute_force(md::system const&, md::array_view<md::vector> forces) override
        {
            for (md::vector& force : forces) {
                force.x += 1;
                force.y += 1;
                force.z += 1;
            }
        }
    };

    md::system system;
    system.add_particle();

    system.add_forcefield(std::make_shared<my_forcefield>());
    system.add_forcefield(std::make_shared<my_forcefield>());

    std::vector<md::vector> forces(system.particle_count());
    system.compute_force(forces);

    CHECK(forces[0].x == 2);
    CHECK(forces[0].y == 2);
    CHECK(forces[0].z == 2);
}
