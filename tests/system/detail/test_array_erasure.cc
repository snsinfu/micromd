#include <algorithm>
#include <memory>

#include <md/system/detail/array_erasure.hpp>

#include <catch.hpp>


TEST_CASE("array_erasure::make - creates an array_erasure with specified size")
{
    using md::detail::array_erasure;

    std::unique_ptr<array_erasure> arr = array_erasure::make<int>(42, -1);

    CHECK(arr->size() == 42);
}

TEST_CASE("array_erasure::recover - recovers the typed implementation")
{
    using md::detail::array_erasure;

    SECTION("int")
    {
        std::unique_ptr<array_erasure> arr = array_erasure::make<int>(42, -1);
        array_erasure::instance<int>& rec = arr->recover<int>();

        CHECK(rec.size() == 42);
        CHECK(rec.data() != nullptr);

        CHECK(std::count(rec.data(), rec.data() + rec.size(), -1) == 42);
    }

    SECTION("double")
    {
        std::unique_ptr<array_erasure> arr = array_erasure::make<double>(42, 1.23);
        array_erasure::instance<double>& rec = arr->recover<double>();

        CHECK(rec.size() == 42);
        CHECK(rec.data() != nullptr);

        CHECK(std::count(rec.data(), rec.data() + rec.size(), 1.23) == 42);
    }
}

TEST_CASE("array_erasure::resize - shortens array keeping values")
{
    using md::detail::array_erasure;

    std::unique_ptr<array_erasure> arr = array_erasure::make<int>(20, -1);
    array_erasure::instance<int>& rec = arr->recover<int>();

    arr->resize(10);

    CHECK(rec.size() == 10);
    CHECK(std::count(rec.data(), rec.data() + rec.size(), -1) == 10);
}

TEST_CASE("array_erasure::resize - extends array with default-valued new elements")
{
    using md::detail::array_erasure;

    std::unique_ptr<array_erasure> arr = array_erasure::make<int>(10, -1);
    array_erasure::instance<int>& rec = arr->recover<int>();

    std::fill_n(rec.data(), rec.size(), 1);
    arr->resize(20);

    CHECK(rec.size() == 20);
    CHECK(std::count(rec.data(), rec.data() + 10, 1) == 10);
    CHECK(std::count(rec.data() + 10, rec.data() + 20, -1) == 10);
}
