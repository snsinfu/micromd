#include <md/system.hpp>
#include <md/system/attribute.hpp>
#include <md/system/particle.hpp>

#include <catch.hpp>


TEST_CASE("particle_ref - provides read access to basic attributes")
{
    md::system system;

    md::basic_particle_data data;
    data.mass = 1.23;
    data.mobility = 4.56;
    data.position = {7, 8, 9};
    data.velocity = {8, 7, 6};
    system.add_particle(data);

    md::particle_ref ref(system, 0);

    CHECK(ref.index == 0);
    CHECK(ref.mass == data.mass);
    CHECK(ref.mobility == data.mobility);
    CHECK((ref.position - data.position).norm() == 0);
    CHECK((ref.velocity - data.velocity).norm() == 0);
}

TEST_CASE("particle_ref - provides write access to basic attributes")
{
    md::system system;
    system.add_particle();

    md::particle_ref ref(system, 0);

    ref.mass = 9.87;
    ref.mobility = 6.54;
    ref.position = {3, 2, 1};
    ref.velocity = {2, 3, 4};

    CHECK(system.view_masses()[0] == 9.87);
    CHECK(system.view_mobilities()[0] == 6.54);
    CHECK((system.view_positions()[0] - md::point{3, 2, 1}).norm() == 0);
    CHECK((system.view_velocities()[0] - md::vector{2, 3, 4}).norm() == 0);
}

TEST_CASE("particle_ref - allows access to any attribute")
{
    md::system system;

    md::attribute_key<md::scalar, struct my_attr_key> my_attr = {};
    system.add_attribute(my_attr);

    md::particle_ref part = system.add_particle();
    part.view(my_attr) = 1.23;

    CHECK(system.view(my_attr)[0] == 1.23);
}

TEST_CASE("particle_iterator - allows access to a range of particles")
{
    md::system system;
    system.add_particle();
    system.add_particle();
    system.add_particle();
    system.add_particle();

    md::particle_iterator it(system, 1);

    md::particle_ref p1 = *it++;
    md::particle_ref p2 = *it++;
    CHECK(it != md::particle_iterator(system, 1));
    CHECK(it == md::particle_iterator(system, 3));

    CHECK(p1.index == 1);
    CHECK(p2.index == 2);
}
