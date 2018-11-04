#include <memory>

#include <md/system.hpp>

#include <catch.hpp>


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

TEST_CASE("system::add_particle - can set basic attributes")
{
    md::system system;

    md::basic_particle_data data;
    data.mass = 1.23;
    data.position = {4, 5, 6};
    data.velocity = {7, 8, 9};
    system.add_particle(data);

    CHECK(system.view(md::mass_attribute)[0] == data.mass);

    CHECK(system.view(md::position_attribute)[0].x == data.position.x);
    CHECK(system.view(md::position_attribute)[0].y == data.position.y);
    CHECK(system.view(md::position_attribute)[0].z == data.position.z);

    CHECK(system.view(md::velocity_attribute)[0].x == data.velocity.x);
    CHECK(system.view(md::velocity_attribute)[0].y == data.velocity.y);
    CHECK(system.view(md::velocity_attribute)[0].z == data.velocity.z);
}

TEST_CASE("system::view_masses - returns mass attribute")
{
    md::system system;

    system.add_particle();
    system.add_particle();
    system.add_particle();

    SECTION("mutable view")
    {
        md::array_view<md::scalar> expected = system.view(md::mass_attribute);
        md::array_view<md::scalar> actual = system.view_masses();

        CHECK(actual.data() == expected.data());
        CHECK(actual.size() == expected.size());
    }

    SECTION("const view")
    {
        md::system const& const_system = system;

        md::array_view<md::scalar const> expected = const_system.view(md::mass_attribute);
        md::array_view<md::scalar const> actual = const_system.view_masses();

        CHECK(actual.data() == expected.data());
        CHECK(actual.size() == expected.size());
    }
}

TEST_CASE("system::view_velocities - returns velocity attribute")
{
    md::system system;

    system.add_particle();
    system.add_particle();
    system.add_particle();

    SECTION("mutable view")
    {
        md::array_view<md::vector> expected = system.view(md::velocity_attribute);
        md::array_view<md::vector> actual = system.view_velocities();

        CHECK(actual.data() == expected.data());
        CHECK(actual.size() == expected.size());
    }

    SECTION("const view")
    {
        md::system const& const_system = system;

        md::array_view<md::vector const> expected = const_system.view(md::velocity_attribute);
        md::array_view<md::vector const> actual = const_system.view_velocities();

        CHECK(actual.data() == expected.data());
        CHECK(actual.size() == expected.size());
    }
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

TEST_CASE("system - has 1-valued mass_attribute by default")
{
    md::system system;
    system.add_particle();

    md::array_view<md::scalar> masses = system.view(md::mass_attribute);

    CHECK(masses.size() == 1);
    CHECK(masses[0] == 1);
}

TEST_CASE("system - has 0-valued position_attribute by default")
{
    md::system system;
    system.add_particle();

    md::array_view<md::point> positions = system.view(md::position_attribute);

    CHECK(positions.size() == 1);
    CHECK(positions[0].x == 0);
    CHECK(positions[0].y == 0);
    CHECK(positions[0].z == 0);
}

TEST_CASE("system - has 0-valued velocity_attribute by default")
{
    md::system system;
    system.add_particle();

    md::array_view<md::vector> velocities = system.view(md::velocity_attribute);

    CHECK(velocities.size() == 1);
    CHECK(velocities[0].x == 0);
    CHECK(velocities[0].y == 0);
    CHECK(velocities[0].z == 0);
}
