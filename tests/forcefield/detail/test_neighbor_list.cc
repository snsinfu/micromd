#include <algorithm>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include <md/forcefield/detail/neighbor_list.hpp>

#include <catch.hpp>


TEST_CASE("neighbor_list - is valid and empty by default")
{
    md::neighbor_list list;

    CHECK(list.begin() == list.end());
}

TEST_CASE("neighbor_list - is valid and empty after updated by zero points")
{
    std::vector<md::point> points;

    md::neighbor_list list;
    list.update(points, 1);

    CHECK(list.begin() == list.end());
}

TEST_CASE("neighbor_list - contains all neighbors in a small example")
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

    std::multiset<std::pair<md::index, md::index>> const expected_pairs = {
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

    md::neighbor_list list;
    list.update(points, cutoff_distance);
    std::multiset<std::pair<md::index, md::index>> pairs(list.begin(), list.end());

    // The list may contain false positives but NO false-negative.
    CHECK(std::includes(
        pairs.begin(),
        pairs.end(),
        expected_pairs.begin(),
        expected_pairs.end()
    ));
}

TEST_CASE("neighbor_list - contains all neighbors on a 3D grid")
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
    md::neighbor_list list;
    list.update(points, cutoff_distance);
    std::multiset<std::pair<md::index, md::index>> pairs(list.begin(), list.end());

    CHECK(std::includes(
        pairs.begin(),
        pairs.end(),
        expected_pairs.begin(),
        expected_pairs.end()
    ));
}

TEST_CASE("neighbor_list - correctly updates after small movement")
{
    int const step_count = 10;
    md::scalar const cutoff_distance = 0.1;
    md::scalar const movement = 0.002;
    md::index const point_count = 50;

    // Generate points
    std::vector<md::point> points;
    std::mt19937 random;
    std::normal_distribution<md::scalar> normal;

    std::generate_n(std::back_inserter(points), point_count, [&] {
        return md::point{normal(random), normal(random), normal(random)};
    });

    md::neighbor_list list;
    list.update(points, cutoff_distance);

    for (int step = 0; step < step_count; step++) {
        // Move points
        for (md::point& pt : points) {
            md::vector delta = {
                normal(random),
                normal(random),
                normal(random)
            };
            delta *= movement / (delta.norm() + movement / 100);

            pt += delta;
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

        // Test on-line update
        list.update(points, cutoff_distance);
        std::multiset<std::pair<md::index, md::index>> pairs(list.begin(), list.end());

        CHECK(std::includes(
            pairs.begin(),
            pairs.end(),
            expected_pairs.begin(),
            expected_pairs.end()
        ));
    }
}

TEST_CASE("neighbor_list - correctly updates after large movement")
{
    int const step_count = 10;
    md::scalar const cutoff_distance = 0.1;
    md::scalar const movement = 0.1;
    md::index const point_count = 50;

    // Generate points
    std::vector<md::point> points;
    std::mt19937 random;
    std::normal_distribution<md::scalar> normal;

    std::generate_n(std::back_inserter(points), point_count, [&] {
        return md::point{normal(random), normal(random), normal(random)};
    });

    md::neighbor_list list;
    list.update(points, cutoff_distance);

    for (int step = 0; step < step_count; step++) {
        // Move points
        for (md::point& pt : points) {
            md::vector delta = {
                normal(random),
                normal(random),
                normal(random)
            };
            delta *= movement / (delta.norm() + movement / 100);

            pt += delta;
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

        // Test on-line update
        list.update(points, cutoff_distance);
        std::multiset<std::pair<md::index, md::index>> pairs(list.begin(), list.end());

        CHECK(std::includes(
            pairs.begin(),
            pairs.end(),
            expected_pairs.begin(),
            expected_pairs.end()
        ));
    }
}

TEST_CASE("neighbor_list - correctly updates after particle count changes")
{
    int const step_count = 5;
    md::index const batch_size = 10;
    md::scalar const cutoff_distance = 0.5;

    md::neighbor_list list;
    std::vector<md::point> points;

    std::mt19937 random;
    std::normal_distribution<md::scalar> normal;

    for (int step = 0; step < step_count; step++) {
        // Generate additional points
        std::generate_n(std::back_inserter(points), batch_size, [&] {
            return md::point{normal(random), normal(random), normal(random)};
        });

        // Brute-force compute neighbor pairs
        std::multiset<std::pair<md::index, md::index>> expected_pairs;

        for (md::index j = 0; j < points.size(); j++) {
            for (md::index i = 0; i < j; i++) {
                if (md::distance(points[i], points[j]) < cutoff_distance) {
                    expected_pairs.emplace(i, j);
                }
            }
        }

        // Test on-line update
        list.update(points, cutoff_distance);
        std::multiset<std::pair<md::index, md::index>> pairs(list.begin(), list.end());

        CHECK(std::includes(
            pairs.begin(),
            pairs.end(),
            expected_pairs.begin(),
            expected_pairs.end()
        ));
    }
}

TEST_CASE("neighbor_list - contains all neighbors in a large cloud", "[.][slow]")
{
    // Generate point cloud
    md::index const point_count = 10000;
    md::scalar const cutoff_distance = 0.1;

    std::mt19937 random;
    std::normal_distribution<md::scalar> normal;

    std::vector<md::point> points;
    std::generate_n(std::back_inserter(points), point_count, [&] {
        return md::point{normal(random), normal(random), normal(random)};
    });

    // Brute-force compute neighbor pairs
    std::multiset<std::pair<md::index, md::index>> expected_pairs;
    for (md::index j = 0; j < points.size(); j++) {
        for (md::index i = 0; i < j; i++) {
            if (md::distance(points[i], points[j]) < cutoff_distance) {
                expected_pairs.emplace(i, j);
            }
        }
    }

    // Test neighbor_list output
    md::neighbor_list list;
    list.update(points, cutoff_distance);
    std::multiset<std::pair<md::index, md::index>> pairs(list.begin(), list.end());

    CHECK(std::includes(
        pairs.begin(),
        pairs.end(),
        expected_pairs.begin(),
        expected_pairs.end()
    ));
}
