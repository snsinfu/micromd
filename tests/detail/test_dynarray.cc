#include <md/detail/dynarray.hpp>

#include "../catch.hpp"


TEST_CASE("dynarray - constructor creates an array")
{
    md::detail::dynarray<int> array{4, 42};

    CHECK(array.view().size() == 4);
    CHECK(array.view()[0] == 42);
    CHECK(array.view()[1] == 42);
    CHECK(array.view()[2] == 42);
    CHECK(array.view()[3] == 42);
}

TEST_CASE("dynarray::resize - extends the contained array with configured initializer")
{
    md::detail::dynarray<int> array{0, 42};

    CHECK(array.view().size() == 0);

    array.resize(2);

    CHECK(array.view().size() == 2);
    CHECK(array.view()[0] == 42);
    CHECK(array.view()[1] == 42);
}

TEST_CASE("dynarray::resize - truncates the contained array")
{
    md::detail::dynarray<int> array{10, 42};

    CHECK(array.view().size() == 10);

    array.resize(2);

    CHECK(array.view().size() == 2);
    CHECK(array.view()[0] == 42);
    CHECK(array.view()[1] == 42);
}

TEST_CASE("dynarray_base - erases dynarray type")
{
    md::detail::dynarray<int> array{0, 42};
    md::detail::dynarray_base& erased = array;

    erased.resize(2);

    auto& recovered = dynamic_cast<md::detail::dynarray<int>&>(erased);

    CHECK(recovered.view().size() == 2);
    CHECK(recovered.view()[0] == 42);
    CHECK(recovered.view()[1] == 42);
}
