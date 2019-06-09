#include <type_traits>

#include <md/misc/math.hpp>

#include <catch.hpp>


TEST_CASE("power - computes correct value")
{
    CHECK(md::power<0>(1.23) == Approx(std::pow(1.23, 0)));
    CHECK(md::power<1>(1.23) == Approx(std::pow(1.23, 1)));
    CHECK(md::power<2>(1.23) == Approx(std::pow(1.23, 2)));
    CHECK(md::power<3>(1.23) == Approx(std::pow(1.23, 3)));
    CHECK(md::power<4>(1.23) == Approx(std::pow(1.23, 4)));
}

TEST_CASE("power_sqrt - computes correct value")
{
    CHECK(md::power_sqrt<0>(1.23) == Approx(std::pow(std::sqrt(1.23), 0)));
    CHECK(md::power_sqrt<1>(1.23) == Approx(std::pow(std::sqrt(1.23), 1)));
    CHECK(md::power_sqrt<2>(1.23) == Approx(std::pow(std::sqrt(1.23), 2)));
    CHECK(md::power_sqrt<3>(1.23) == Approx(std::pow(std::sqrt(1.23), 3)));
    CHECK(md::power_sqrt<4>(1.23) == Approx(std::pow(std::sqrt(1.23), 4)));
}
