#include <md/basic_types.hpp>
#include <md/system/attribute.hpp>
#include <md/system/detail/attribute_table.hpp>

#include <catch.hpp>


TEST_CASE("attribute_table - is default constructible and empty")
{
    md::detail::attribute_table table;

    CHECK(table.size() == 0);
}

TEST_CASE("attribute_table::resize - changes the size of table")
{
    md::detail::attribute_table table;

    table.resize(10);
    CHECK(table.size() == 10);

    table.resize(20);
    CHECK(table.size() == 20);

    table.resize(5);
    CHECK(table.size() == 5);
}

TEST_CASE("attribute_table::require - allocates a column")
{
    md::detail::attribute_table table;

    md::attribute_key<int, struct test> test_attribute = [](test*) -> int {
        return 42;
    };

    table.require(test_attribute);
    table.resize(4);

    md::array_view<int> values = table.view(test_attribute);

    CHECK(values.size() == 4);
    CHECK(values[0] == 42);
    CHECK(values[1] == 42);
    CHECK(values[2] == 42);
    CHECK(values[3] == 42);
}

TEST_CASE("attribute_table::require - is idempotent")
{
    md::detail::attribute_table table;

    md::attribute_key<int, struct test> test_attribute = [](test*) -> int {
        return 42;
    };

    table.require(test_attribute);
    table.resize(4);

    md::array_view<int> values = table.view(test_attribute);

    CHECK(values.size() == 4);
    CHECK(values[0] == 42);
    CHECK(values[1] == 42);
    CHECK(values[2] == 42);
    CHECK(values[3] == 42);

    table.require(test_attribute);

    CHECK(values.size() == 4);
    CHECK(values[0] == 42);
    CHECK(values[1] == 42);
    CHECK(values[2] == 42);
    CHECK(values[3] == 42);
}

TEST_CASE("attribute_table - is copyable")
{
    md::detail::attribute_table table;

    md::attribute_key<int, struct test> test_attribute = [](test*) -> int {
        return 42;
    };

    table.require(test_attribute);
    table.resize(4);

    SECTION("copy construction")
    {
        md::detail::attribute_table copy = table;
        CHECK(copy.size() == 4);

        md::array_view<int> values = copy.view(test_attribute);
        CHECK(values.size() == 4);
        CHECK(values[0] == 42);
        CHECK(values[1] == 42);
        CHECK(values[2] == 42);
        CHECK(values[3] == 42);
    }

    SECTION("copy assignment")
    {
        md::detail::attribute_table copy;
        copy = table;
        CHECK(copy.size() == 4);

        md::array_view<int> values = copy.view(test_attribute);
        CHECK(values.size() == 4);
        CHECK(values[0] == 42);
        CHECK(values[1] == 42);
        CHECK(values[2] == 42);
        CHECK(values[3] == 42);
    }
}
