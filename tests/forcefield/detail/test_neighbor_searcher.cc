#include <algorithm>
#include <cmath>
#include <iterator>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include <md/forcefield/detail/neighbor_searcher.hpp>

#include <catch.hpp>


TEST_CASE("neighbor_searcher - outputs nothing by default")
{
    std::vector<std::pair<md::index, md::index>> pairs;
    md::neighbor_searcher searcher(1, md::linear_hash{});
    searcher.search(std::back_inserter(pairs));

    CHECK(pairs.empty());
}

TEST_CASE("neighbor_searcher - outputs nothing for empty points")
{
    std::vector<md::point> points;
    std::vector<std::pair<md::index, md::index>> pairs;

    md::neighbor_searcher searcher(1, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::back_inserter(pairs));

    CHECK(pairs.empty());
}

TEST_CASE("neighbor_searcher - finds neighbors on a line")
{
    md::scalar const cutoff_distance = 0.1;
    md::vector const step_vector = {0.06571881, 0.00727112, 0.02298191};

    // Generate points on a line
    std::vector<md::point> points;
    std::multiset<std::pair<md::index, md::index>> expected_pairs;

    points.push_back({-1, -1, -1});
    for (md::index i = 1; i < 50; i++) {
        points.push_back(points.front() + md::scalar(i) * step_vector);
        expected_pairs.emplace(i - 1, i);
    }
    CHECK(expected_pairs.size() == points.size() - 1);

    // Test
    std::multiset<std::pair<md::index, md::index>> pairs;
    md::neighbor_searcher searcher(cutoff_distance, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::inserter(pairs, pairs.end()));

    CHECK(pairs.size() == expected_pairs.size());
    CHECK(pairs == expected_pairs);
}

TEST_CASE("neighbor_searcher - finds neighbors on a 2D grid")
{
    md::scalar const cutoff_distance = 0.1;
    md::scalar const spacing = 0.07;

    // Generate 2D grid points
    std::vector<md::point> points;

    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            points.push_back({
                spacing * x,
                spacing * y,
                0
            });
        }
    }

    // Brute-force compute neighbor pairs
    std::multiset<std::pair<md::index, md::index>> expected_pairs;

    for (md::index j = 0; j < points.size(); j++) {
        for (md::index i = 0; i < j; i++) {
            if (md::distance(points[i], points[j]) < cutoff_distance) {
                expected_pairs.emplace(i, j);
            }
        }
    }

    // Test
    std::multiset<std::pair<md::index, md::index>> pairs;
    md::neighbor_searcher searcher(cutoff_distance, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::inserter(pairs, pairs.end()));

    CHECK(pairs.size() == expected_pairs.size());
    CHECK(pairs == expected_pairs);
}

TEST_CASE("neighbor_searcher - finds neighbors on a 3D grid")
{
    md::scalar const cutoff_distance = 0.1;
    md::scalar const spacing = 0.07;

    // Generate a 3D grid points
    std::vector<md::point> points;

    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            for (int z = -2; z <= 2; z++) {
                points.push_back({
                    spacing * x,
                    spacing * y,
                    spacing * z
                });
            }
        }
    }

    // Brute-force compute neighbor pairs
    std::multiset<std::pair<md::index, md::index>> expected_pairs;

    for (md::index j = 0; j < points.size(); j++) {
        for (md::index i = 0; i < j; i++) {
            if (md::distance(points[i], points[j]) < cutoff_distance) {
                expected_pairs.emplace(i, j);
            }
        }
    }

    // Test
    std::multiset<std::pair<md::index, md::index>> pairs;
    md::neighbor_searcher searcher(cutoff_distance, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::inserter(pairs, pairs.end()));

    CHECK(pairs.size() == expected_pairs.size());
    CHECK(pairs == expected_pairs);
}

TEST_CASE("neighbor_searcher - finds neighbors in a random walk path")
{
    md::scalar const cutoff_distance = 0.1;
    md::scalar const spacing = 0.07;

    // Generate a random walk path
    std::vector<md::point> points;

    std::mt19937_64 random;
    std::normal_distribution<md::scalar> normal;

    points.push_back({0, 0, 0});
    for (int i = 1; i < 50; i++) {
        md::vector step = {normal(random), normal(random), normal(random)};
        step *= spacing / step.norm();
        points.push_back(points.back() + step);
    }

    // Brute-force computed neighbor pairs
    std::multiset<std::pair<md::index, md::index>> expected_pairs;

    for (md::index j = 0; j < points.size(); j++) {
        for (md::index i = 0; i < j; i++) {
            if (md::distance(points[i], points[j]) < cutoff_distance) {
                expected_pairs.emplace(i, j);
            }
        }
    }

    // Test
    std::multiset<std::pair<md::index, md::index>> pairs;
    md::neighbor_searcher searcher(cutoff_distance, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::inserter(pairs, pairs.end()));

    CHECK(pairs.size() == expected_pairs.size());
    CHECK(pairs == expected_pairs);
}

TEST_CASE("neighbor_searcher - finds no false positives in a sparse system")
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
    md::neighbor_searcher searcher(cutoff_distance, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::inserter(actual, actual.end()));

    CHECK(actual.empty());
}

TEST_CASE("neighbor_searcher - finds correct neighbors in a small example")
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

    md::neighbor_searcher searcher(cutoff_distance, hash);
    searcher.set_points(points);
    searcher.search(std::inserter(actual, actual.end()));

    CHECK(actual == expected);
}

TEST_CASE("neighbor_searcher - finds correct neighbors in a large cloud", "[.][slow]")
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
    md::neighbor_searcher searcher(cutoff_distance, md::linear_hash{});
    searcher.set_points(points);
    searcher.search(std::inserter(actual, actual.end()));

    CHECK(actual == expected);
}

TEST_CASE("neighbor_searcher - finds neighbors of query point in a point cloud")
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

    md::linear_hash hash;
    hash.x_coeff = 220598453;
    hash.y_coeff = 126390227;
    hash.z_coeff = 435046583;
    hash.modulus = 11;

    md::neighbor_searcher searcher(cutoff_distance, hash);
    searcher.set_points(points);

    SECTION("example 1")
    {
        md::point const query_point = {-0.3, 0.6, 0.1};
        std::multiset<md::index> const expected = {4, 7};
        std::multiset<md::index> actual;

        searcher.query(query_point, std::inserter(actual, actual.end()));

        CHECK(actual == expected);
    }

    SECTION("example 2")
    {
        md::point const query_point = {1.0, 0.0, 0.0};
        std::multiset<md::index> const expected = {0, 5, 6, 9};
        std::multiset<md::index> actual;

        searcher.query(query_point, std::inserter(actual, actual.end()));

        CHECK(actual == expected);
    }

    SECTION("example 3: no neighbor")
    {
        md::point const query_point = {-1.0, 1.0, 1.0};
        std::multiset<md::index> actual;

        searcher.query(query_point, std::inserter(actual, actual.end()));

        CHECK(actual.empty());
    }
}
