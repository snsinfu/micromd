#include <type_traits>

#include "../include/md.hpp"
#include "catch.hpp"


TEST_CASE("linear_hash - hash type is unsigned")
{
    CHECK(std::is_unsigned<md::linear_hash::hash_t>::value);
}

TEST_CASE("linear_hash - is default constructible and modulus is nonzero")
{
    md::linear_hash hash;

    CHECK(hash.modulus > 0);
}

TEST_CASE("linear_hash::operator() - trivial zero points")
{
    md::linear_hash hash;

    CHECK(hash(0, 0, 0) == 0);
    CHECK(hash(hash.modulus, hash.modulus, hash.modulus) == 0);
}

TEST_CASE("linear_hash::operator() - computes expected hash value")
{
    md::linear_hash hash;

    hash.x_coeff = 2;
    hash.y_coeff = 3;
    hash.z_coeff = 5;
    hash.modulus = 11;

    md::linear_hash const& const_hash = hash;

    CHECK(const_hash(1, 0, 0) == 2);
    CHECK(const_hash(0, 1, 0) == 3);
    CHECK(const_hash(0, 0, 1) == 5);
    CHECK(const_hash(1, 1, 1) == 10);

    CHECK(const_hash(1000, 10, 100) == 0);
    CHECK(const_hash(1000, 100, 10) == 7);
    CHECK(const_hash(10, 1000, 100) == 0);
    CHECK(const_hash(100, 1000, 10) == 5);
    CHECK(const_hash(10, 100, 1000) == 7);
    CHECK(const_hash(100, 10, 1000) == 5);
}

TEST_CASE("linear_hash::operator() - avoids 32-bit wraparound")
{
    md::linear_hash hash;

    hash.x_coeff = 3416542937;
    hash.y_coeff = 3481543853;
    hash.z_coeff = 3912176573;
    hash.modulus = 17;

    CHECK(hash(1234, 5678, 9012) == 14);
    CHECK(hash(3456, 7890, 1234) == 4);
    CHECK(hash(5678, 9012, 3456) == 0);
}
