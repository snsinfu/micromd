#include <md/system/attribute.hpp>

#include <catch.hpp>


TEST_CASE("attribute_key - is a function pointer type")
{
    md::attribute_key<int, struct test> key = +[](test*) -> int {
        return 0;
    };

    CHECK(key);
}

TEST_CASE("default_value - returns a default constructed value for a null key")
{
    md::attribute_key<int, struct test> key = {};

    CHECK(md::default_value(key) == 0);
}

TEST_CASE("default_value - returns a function result for a non-null key")
{
    md::attribute_key<int, struct test> key = [](test*) -> int {
        return 42;
    };

    CHECK(md::default_value(key) == 42);
}
