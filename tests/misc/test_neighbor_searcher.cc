#include <algorithm>
#include <cmath>
#include <iterator>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include <md/misc/box.hpp>
#include <md/misc/neighbor_searcher.hpp>

#include <catch.hpp>


namespace
{
    // Returns true if the multiset contains the same values of the other set.
    template<typename T>
    bool equal(std::multiset<T> const& set1, std::set<T> const& set2)
    {
        return std::equal(set1.begin(), set1.end(), set2.begin(), set2.end());
    }
}


TEST_CASE("neighbor_searcher - is instantiable with standard boxes")
{
    using open_searcher = md::neighbor_searcher<md::open_box>;
    using periodic_searcher = md::neighbor_searcher<md::periodic_box>;
    using xy_periodic_searcher = md::neighbor_searcher<md::xy_periodic_box>;
    CHECK(sizeof(open_searcher) > 0);
    CHECK(sizeof(periodic_searcher) > 0);
    CHECK(sizeof(xy_periodic_searcher) > 0);
}

TEST_CASE("neighbor_searcher - outputs nothing by default")
{
    md::neighbor_searcher<md::open_box> searcher{md::open_box{}, 1};
    std::vector<std::pair<md::index, md::index>> pairs;
    searcher.search(std::back_inserter(pairs));
    CHECK(pairs.empty());
}

TEST_CASE("neighbor_searcher - accepts empty point cloud")
{
    md::neighbor_searcher<md::open_box> searcher{md::open_box{}, 1};
    std::vector<std::pair<md::index, md::index>> pairs;
    std::vector<md::point> points;
    searcher.set_points(points);
    searcher.search(std::back_inserter(pairs));
    CHECK(pairs.empty());
}

TEST_CASE("neighbor_searcher - finds correct neighbors in point cloud")
{
    md::index const point_count = 100;
    md::scalar const neighbor_distance = 0.3;

    std::mt19937 random;
    std::normal_distribution<md::scalar> coord;
    std::vector<md::point> points;
    std::generate_n(std::back_inserter(points), point_count, [&] {
        return md::point{coord(random), coord(random), coord(random)};
    });

    struct actual_expect
    {
        // Computed neighbor pairs. This is multiset, not set, so that errornous
        // duplicate pairs are catched.
        std::multiset<std::pair<md::index, md::index>> actual;

        // Ground truth.
        std::set<std::pair<md::index, md::index>> expect;
    };

    auto compute_actual_expect = [&](auto box) {
        actual_expect result;

        // Ground truth by brute-force method.
        for (md::index j = 0; j < points.size(); j++) {
            for (md::index i = 0; i < j; i++) {
                md::vector const disp = box.shortest_displacement(points[i], points[j]);
                md::scalar const dist = disp.norm();
                if (dist < neighbor_distance) {
                    result.expect.emplace(i, j);
                }
            }
        }

        // Actual output from neighbor_searcher.
        using box_type = decltype(box);
        md::neighbor_searcher<box_type> searcher{box, neighbor_distance};
        searcher.set_points(points);
        searcher.search(std::inserter(result.actual, result.actual.end()));

        return result;
    };

    SECTION("open_box")
    {
        md::open_box box;
        box.particle_count = point_count;
        auto const result = compute_actual_expect(box);
        CHECK(equal(result.actual, result.expect));
    }

    SECTION("periodic_box")
    {
        md::periodic_box box;
        box.x_period = 0.9;
        box.y_period = 1.0;
        box.z_period = 1.1;
        auto const result = compute_actual_expect(box);
        CHECK(equal(result.actual, result.expect));
    }

    SECTION("xy_periodic_box")
    {
        md::xy_periodic_box box;
        box.x_period = 0.9;
        box.y_period = 1.1;
        box.z_span = 1.0;
        box.particle_count = point_count;
        auto const result = compute_actual_expect(box);
        CHECK(equal(result.actual, result.expect));
    }
}

TEST_CASE("neighbor_searcher - finds correct neighbors of query point")
{
    md::index const point_count = 100;
    md::scalar const neighbor_distance = 0.3;

    std::mt19937 random;
    std::normal_distribution<md::scalar> coord;
    std::vector<md::point> points;
    std::generate_n(std::back_inserter(points), point_count, [&] {
        return md::point{coord(random), coord(random), coord(random)};
    });

    struct actual_expect
    {
        // Computed neighbor pairs. This is multiset, not set, so that errornous
        // duplicate pairs are catched.
        std::multiset<md::index> actual;

        // Ground truth.
        std::set<md::index> expect;
    };

    auto compute_actual_expect = [&](auto box, md::point query) {
        actual_expect result;

        // Ground truth by brute-force method.
        for (md::index i = 0; i < points.size(); i++) {
            md::vector const disp = box.shortest_displacement(points[i], query);
            md::scalar const dist = disp.norm();
            if (dist < neighbor_distance) {
                result.expect.insert(i);
            }
        }

        // Actual output from neighbor_searcher.
        using box_type = decltype(box);
        md::neighbor_searcher<box_type> searcher{box, neighbor_distance};
        searcher.set_points(points);
        searcher.query(query, std::inserter(result.actual, result.actual.end()));

        return result;
    };

    SECTION("open_box")
    {
        md::open_box box;
        box.particle_count = point_count;
        auto const result = compute_actual_expect(box, {-0.1, 0.0, 0.1});
        CHECK(equal(result.actual, result.expect));
    }

    SECTION("periodic_box")
    {
        md::periodic_box box;
        box.x_period = 0.9;
        box.y_period = 1.0;
        box.z_period = 1.1;
        auto const result = compute_actual_expect(box, {-0.1, 0.0, 0.1});
        CHECK(equal(result.actual, result.expect));
    }

    SECTION("xy_periodic_box")
    {
        md::xy_periodic_box box;
        box.x_period = 0.9;
        box.y_period = 1.1;
        box.z_span = 1.0;
        box.particle_count = point_count;
        auto const result = compute_actual_expect(box, {-0.1, 0.0, 0.1});
        CHECK(equal(result.actual, result.expect));
    }
}
