#include <md/basic_types.hpp>
#include <md/misc/box.hpp>

#include <catch.hpp>


TEST_CASE("open_box - has reasonable defaults for hint parameters")
{
    md::open_box box;
    CHECK(box.particle_count > 0);
}

TEST_CASE("xy_periodic_box - has reasonable defaults for hint parameters")
{
    md::xy_periodic_box box;
    CHECK(box.z_span > 0);
    CHECK(box.particle_count > 0);
}

TEST_CASE("open_box - computes correct displacement")
{
    md::open_box box;

    md::point const p1 = {1, 2, 3};
    md::point const p2 = {4, 5, 6};
    md::vector const expected = p1 - p2;
    md::vector const actual = box.shortest_displacement(p1, p2);
    CHECK(actual.x == Approx(expected.x));
    CHECK(actual.y == Approx(expected.y));
    CHECK(actual.z == Approx(expected.z));
}

TEST_CASE("periodic_box - computes correct displacement")
{
    md::periodic_box box;
    box.x_period = 1;
    box.y_period = 2;
    box.z_period = 3;

    SECTION("internal")
    {
        md::point const p1 = {0.4, 0.5, 0.6};
        md::point const p2 = {0.6, 0.8, 1.0};
        md::vector const expected = {-0.2, -0.3, -0.4};
        md::vector const actual = box.shortest_displacement(p1, p2);
        CHECK(actual.x == Approx(expected.x));
        CHECK(actual.y == Approx(expected.y));
        CHECK(actual.z == Approx(expected.z));
    }

    SECTION("edge to edge")
    {
        md::point const p1 = {0.1, 0.2, 0.3};
        md::point const p2 = {0.9, 1.9, 2.9};
        md::vector const expected = {0.2, 0.3, 0.4};
        md::vector const actual = box.shortest_displacement(p1, p2);
        CHECK(actual.x == Approx(expected.x));
        CHECK(actual.y == Approx(expected.y));
        CHECK(actual.z == Approx(expected.z));
    }
}

TEST_CASE("xy_periodic_box - computes correct displacement")
{
    md::xy_periodic_box box;
    box.x_period = 1;
    box.y_period = 2;

    SECTION("internal")
    {
        md::point const p1 = {0.4, 0.5, 0.6};
        md::point const p2 = {0.6, 0.8, 1.0};
        md::vector const expected = {-0.2, -0.3, -0.4};
        md::vector const actual = box.shortest_displacement(p1, p2);
        CHECK(actual.x == Approx(expected.x));
        CHECK(actual.y == Approx(expected.y));
        CHECK(actual.z == Approx(expected.z));
    }

    SECTION("edge to edge")
    {
        md::point const p1 = {0.1, 0.2, 0.3};
        md::point const p2 = {0.9, 1.9, 2.9};
        md::vector const expected = {0.2, 0.3, -2.6};
        md::vector const actual = box.shortest_displacement(p1, p2);
        CHECK(actual.x == Approx(expected.x));
        CHECK(actual.y == Approx(expected.y));
        CHECK(actual.z == Approx(expected.z));
    }
}
