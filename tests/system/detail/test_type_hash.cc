#include <algorithm>
#include <memory>

#include <md/system/detail/type_hash.hpp>

#include <catch.hpp>


TEST_CASE("type_hash::hash_v - is unique")
{
    using md::detail::type_hash;

    std::size_t const foo = type_hash::hash<struct foo_tag>::value;
    std::size_t const bar = type_hash::hash<struct bar_tag>::value;
    std::size_t const baz = type_hash::hash<struct baz_tag>::value;

    CHECK(foo != bar);
    CHECK(bar != baz);
    CHECK(baz != foo);
}

TEST_CASE("type_hash::hash_v - is less than size()")
{
    using md::detail::type_hash;

    std::size_t const foo = type_hash::hash<struct foo_tag>::value;
    std::size_t const bar = type_hash::hash<struct bar_tag>::value;
    std::size_t const baz = type_hash::hash<struct baz_tag>::value;

    CHECK(foo < type_hash::size());
    CHECK(bar < type_hash::size());
    CHECK(baz < type_hash::size());
}
