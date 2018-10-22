#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "../include/md.hpp"
#include "catch.hpp"


TEST_CASE("hashing_neighbor_searcher - outputs nothing by default")
{
    std::vector<std::pair<md::index, md::index>> pairs;
    md::hashing_neighbor_searcher searcher(1, md::linear_hash{});
    searcher.search(std::back_inserter(pairs));

    CHECK(pairs.empty());
}

TEST_CASE("hashing_neighbor_searcher - outputs nothing for empty points")
{
    std::vector<md::point> points;
    std::vector<std::pair<md::index, md::index>> pairs;

    md::hashing_neighbor_searcher searcher(1, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::back_inserter(pairs));

    CHECK(pairs.empty());
}

TEST_CASE("hashing_neighbor_searcher - finds correct neighbors")
{
    md::scalar const cutoff_distance = 1.0;

    std::vector<md::point> const points = {
        { 0.6, -0.5,  0.4},
        { 0.9, -1.0, -0.2},
        {-0.5, -1.0,  0.0},
        {-0.3, -0.5,  0.1},
        {-0.7,  0.1,  0.0},
        { 0.7,  0.0, -0.2},
        { 0.5, -0.5, -0.5},
        {-0.4,  0.8, -0.7},
        {-0.5, -0.6,  0.4},
        { 0.9,  0.5, -0.8},
    };

    std::multiset<std::pair<md::index, md::index>> const expected = {
        {0, 1},
        {0, 3},
        {0, 5},
        {0, 6},
        {1, 6},
        {2, 3},
        {2, 8},
        {3, 4},
        {3, 6},
        {3, 8},
        {4, 8},
        {5, 6},
        {5, 9},
    };

    std::multiset<std::pair<md::index, md::index>> actual;

    md::linear_hash hash;
    hash.x_coeff = 220598453;
    hash.y_coeff = 126390227;
    hash.z_coeff = 435046583;
    hash.modulus = 11;

    md::hashing_neighbor_searcher searcher(cutoff_distance, hash);
    searcher.set_points(points);
    searcher.search(std::inserter(actual, actual.end()));

    CHECK(actual == expected);
}

TEST_CASE("hashing_neighbor_searcher - detects no false positive")
{
    md::scalar const cutoff_distance = 0.1;

    std::vector<md::point> points;

    for (int i = 0; i < 1000; i++) {
        int const x = i % 10;
        int const y = i / 10 % 10;
        int const z = i / 10 / 10 % 10;

        points.push_back({0.11 * x, 0.11 * y, 0.11 * z});
    }

    std::set<std::pair<md::index, md::index>> actual;
    md::hashing_neighbor_searcher searcher(cutoff_distance, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::inserter(actual, actual.end()));

    CHECK(actual.empty());
}

TEST_CASE("hashing_neighbor_searcher - normal distribution", "[.][slow]")
{
    std::mt19937 random;
    std::normal_distribution<md::scalar> normal;

    md::index const point_count = 10000;
    md::scalar const cutoff_distance = 0.1;

    std::vector<md::point> points;
    std::generate_n(std::back_inserter(points), point_count, [&] {
        return md::point{normal(random), normal(random), normal(random)};
    });

    std::multiset<std::pair<md::index, md::index>> expected;
    for (md::index j = 0; j < points.size(); j++) {
        for (md::index i = 0; i < j; i++) {
            if (md::distance(points[i], points[j]) < cutoff_distance) {
                expected.emplace(i, j);
            }
        }
    }

    std::multiset<std::pair<md::index, md::index>> actual;
    md::hashing_neighbor_searcher searcher(cutoff_distance, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::inserter(actual, actual.end()));

    CHECK(actual == expected);
}

TEST_CASE("detail::compute_variance - computes correct variance")
{
    std::vector<md::point> const points = {
        {-0.3,  0.4,  -9.4},
        { 0.9, -0.6,   0.8},
        {-1.7,  0.9,   5.3},
        { 0.0,  3.5,  -2.7},
        {-2.1,  2.0,  -3.1},
        { 3.0, -2.5, -12.8},
        {-1.8, -1.6,   5.1},
        {-1.3, -1.4,   1.7},
        {-0.2,  6.3,  -1.9},
        {-0.6, -6.9,  -5.3}
    };

    md::vector const var = md::detail::compute_variance(points);

    CHECK(var.x == Approx(2.0849));
    CHECK(var.y == Approx(11.5649));
    CHECK(var.z == Approx(30.8701));
}

TEST_CASE("detail::determine_hash - squashes lowest-variance axis")
{
    std::vector<md::point> const points = {
        {-0.3,  0.4,  -9.4},
        { 0.9, -0.6,   0.8},
        {-1.7,  0.9,   5.3},
        { 0.0,  3.5,  -2.7},
        {-2.1,  2.0,  -3.1},
        { 3.0, -2.5, -12.8},
        {-1.8, -1.6,   5.1},
        {-1.3, -1.4,   1.7},
        {-0.2,  6.3,  -1.9},
        {-0.6, -6.9,  -5.3}
    };

    md::linear_hash const hash = md::detail::determine_hash(points);

    // x axis has the lowest variance
    CHECK(hash.x_coeff == 0);
    CHECK(hash.y_coeff != 0);
    CHECK(hash.z_coeff != 0);
    CHECK(hash.modulus != 0);
}
