#include <iterator>

#include <md/basic_types.hpp>
#include <md/misc/index_range.hpp>

#include <catch.hpp>


TEST_CASE("index_iterator - is default constructible and pointing to zero")
{
    md::index_iterator it;
    CHECK(*it == 0);
}

TEST_CASE("index_iterator - is constructible from a number")
{
    md::index_iterator it{10};
    CHECK(*it == 10);
}

TEST_CASE("index_iterator - is equality comparable")
{
    md::index_iterator it_1{10};
    md::index_iterator it_2{20};
    CHECK(it_1 == it_1);
    CHECK(it_1 != it_2);
    CHECK(it_2 == it_2);
}

TEST_CASE("index_iterator - supports incrementing")
{
    SECTION("post-increment")
    {
        md::index_iterator it{10};
        CHECK(*it++ == 10);
        CHECK(*it == 11);
    }

    SECTION("pre-increment")
    {
        md::index_iterator it{10};
        CHECK(*++it == 11);
        CHECK(*it == 11);
    }
}

TEST_CASE("index_range - is default constructible")
{
    md::index_range indices;
    CHECK(indices.size() == 0);
    CHECK(indices.begin() == indices.end());
}

TEST_CASE("index_range - is constructible from a begin-end pair")
{
    md::index_range indices{5, 15};
    CHECK(indices.size() == 10);
    CHECK(*indices.begin() == 5);

    indices = {8, 11};
    CHECK(indices.size() == 3);
    CHECK(*indices.begin() == 8);
}

TEST_CASE("index_range - is constructible from a size number")
{
    md::index_range indices{10};
    CHECK(indices.size() == 10);
    CHECK(*indices.begin() == 0);
}

TEST_CASE("index_range - is a range")
{
    md::index_range indices{3};
    for (md::index i : indices) {
        CHECK(i < 3);
    }
}

TEST_CASE("index_range - supports random access")
{
    md::index_range indices{5, 10};
    CHECK(indices[0] == 5);
    CHECK(indices[3] == 8);
}
