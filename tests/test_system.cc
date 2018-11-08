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
    data.mobility = 3.21;
    data.position = {4, 5, 6};
    data.velocity = {7, 8, 9};
    system.add_particle(data);

    CHECK(system.view_masses()[0] == data.mass);

    CHECK(system.view_mobilities()[0] == data.mobility);

    CHECK(system.view_positions()[0].x == data.position.x);
    CHECK(system.view_positions()[0].y == data.position.y);
    CHECK(system.view_positions()[0].z == data.position.z);

    CHECK(system.view_velocities()[0].x == data.velocity.x);
    CHECK(system.view_velocities()[0].y == data.velocity.y);
    CHECK(system.view_velocities()[0].z == data.velocity.z);
}

TEST_CASE("system::add_particle - returns a particle reference")
{
    md::system system;

    md::particle_ref part = system.add_particle();
    part.mass = 1.23;
    part.mobility = 4.56;
    part.position = {7, 8, 9};
    part.velocity = {8, 7, 6};

    CHECK(system.view_masses()[0] == 1.23);

    CHECK(system.view_mobilities()[0] == 4.56);

    CHECK(system.view_positions()[0].x == 7);
    CHECK(system.view_positions()[0].y == 8);
    CHECK(system.view_positions()[0].z == 9);

    CHECK(system.view_velocities()[0].x == 8);
    CHECK(system.view_velocities()[0].y == 7);
    CHECK(system.view_velocities()[0].z == 6);
}

TEST_CASE("system::particles - returns a range of particle references")
{
    md::system system;

    system.add_particle().mass = 42;
    system.add_particle().mass = 42;
    system.add_particle().mass = 42;

    for (md::particle_ref part : system.particles()) {
        CHECK(part.mass == 42);
    }
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

TEST_CASE("system::view_mobilities - returns mobility attribute")
{
    md::system system;

    system.add_particle();
    system.add_particle();
    system.add_particle();

    SECTION("mutable view")
    {
        md::array_view<md::scalar> expected = system.view(md::mobility_attribute);
        md::array_view<md::scalar> actual = system.view_mobilities();

        CHECK(actual.data() == expected.data());
        CHECK(actual.size() == expected.size());
    }

    SECTION("const view")
    {
        md::system const& const_system = system;

        md::array_view<md::scalar const> expected = const_system.view(md::mobility_attribute);
        md::array_view<md::scalar const> actual = const_system.view_mobilities();

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

TEST_CASE("system::add_forcefield - accepts a forcefield and returns a shared_ptr")
{
    class my_forcefield : public md::forcefield
    {
    public:
        md::scalar energy;

        explicit my_forcefield(md::scalar energy)
            : energy{energy}
        {
        }

        md::scalar compute_energy(md::system const&) override
        {
            return energy;
        }

        void compute_force(md::system const&, md::array_view<md::vector>) override
        {
        }
    };

    md::system system;

    // system wraps a copy of the passed-in ff as a shared_ptr
    std::shared_ptr<my_forcefield> ff = system.add_forcefield(my_forcefield{1.23});

    CHECK(ff->energy == 1.23);

    // The forcefield object is shared
    ff->energy = 3.21;

    CHECK(system.compute_potential_energy() == 3.21);
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

TEST_CASE("system::compute_kinetic_energy - returns the sum of particle kinetic energy")
{
    md::system system;

    md::basic_particle_data part1;
    part1.mass = 1.2;
    part1.velocity = {3.4, 5.6, 7.8};

    md::basic_particle_data part2;
    part2.mass = 9.0;
    part2.velocity = {1.2, 3.4, 5.6};

    system.add_particle(part1);
    system.add_particle(part2);

    md::scalar const kin1 = part1.mass * part1.velocity.squared_norm() / 2;
    md::scalar const kin2 = part2.mass * part2.velocity.squared_norm() / 2;
    md::scalar const expected_energy = kin1 + kin2;

    CHECK(system.compute_kinetic_energy() == Approx(expected_energy));
}

TEST_CASE("system::compute_energy - returns kinetic + potential energy")
{
    class forcefield : public md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const&) override
        {
            return 2.35;
        }

        void compute_force(md::system const&, md::array_view<md::vector>) override
        {
        }
    };

    md::system system;

    md::basic_particle_data part1;
    part1.mass = 1.2;
    part1.velocity = {3.4, 5.6, 7.8};

    md::basic_particle_data part2;
    part2.mass = 9.0;
    part2.velocity = {1.2, 3.4, 5.6};

    system.add_particle(part1);
    system.add_particle(part2);

    system.add_forcefield(std::make_shared<forcefield>());

    md::scalar const kin1 = part1.mass * part1.velocity.squared_norm() / 2;
    md::scalar const kin2 = part2.mass * part2.velocity.squared_norm() / 2;
    md::scalar const pot = 2.35;
    md::scalar const expected_energy = kin1 + kin2 + pot;

    CHECK(system.compute_energy() == Approx(expected_energy));
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
