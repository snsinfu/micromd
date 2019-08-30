#include <algorithm>
#include <iterator>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include <md/forcefield/detail/neighbor_list_v2.hpp>
#include <md/misc/box.hpp>

#include <catch.hpp>


TEST_CASE("neighbor_list_v2 - is valid and empty by default")
{
    md::neighbor_list_v2<md::open_box> list;

    CHECK(list.begin() == list.end());
}

TEST_CASE("neighbor_list_v2 - is valid and empty after updated by zero points")
{
    std::vector<md::point> points;

    md::neighbor_list_v2<md::open_box> list;
    list.update(points, 1, md::open_box{});

    CHECK(list.begin() == list.end());
}

TEST_CASE("neighbor_list_v2 - finds correct neighbor pairs")
{
    md::scalar const cutoff_distance = 0.1;
    md::index const point_count = 1000;

    std::vector<md::point> points;
    std::mt19937 random;
    std::generate_n(std::back_inserter(points), point_count, [&] {
        std::uniform_real_distribution<md::scalar> coord;
        return md::point{coord(random), coord(random), coord(random)};
    });

    auto test_on_box = [&](auto box) {
        // Expect: brute-force.
        std::set<std::pair<md::index, md::index>> expect;

        for (md::index i = 0; i < points.size(); i++) {
            for (md::index j = i + 1; j < points.size(); j++) {
                if (md::distance(points[i], points[j]) < cutoff_distance) {
                    expect.emplace(i, j);
                }
            }
        }

        // Actual: neighbor_list.
        std::set<std::pair<md::index, md::index>> actual;

        using box_type = decltype(box);
        md::neighbor_list_v2<box_type> list;
        list.update(points, cutoff_distance, box);
        for (auto pair : list) {
            actual.insert(pair);
        }

        // The list may contain false positives (for efficient list reuse), but
        // it should never contain any false-negatives. So test actual ⊇ expect.
        CHECK(std::includes(
            actual.begin(), actual.end(), expect.begin(), expect.end()
        ));
    };

    SECTION("in open_box")
    {
        md::open_box box;
        test_on_box(box);
    }

    SECTION("in periodic_box")
    {
        md::periodic_box box;
        box.x_period = 0.9;
        box.y_period = 1.0;
        box.z_period = 1.1;
        test_on_box(box);
    }

    SECTION("in xy_periodic_box")
    {
        md::xy_periodic_box box;
        box.x_period = 0.9;
        box.y_period = 1.0;
        test_on_box(box);
    }
}
