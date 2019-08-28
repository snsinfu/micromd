#include <algorithm>
#include <cmath>
#include <iterator>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include <md/misc/box.hpp>
#include <md/misc/neighbor_searcher_v2.hpp>

#include <catch.hpp>


TEST_CASE("neighbor_searcher_v2 - is instantiable with standard boxes")
{
    using open_searcher = md::neighbor_searcher_v2<md::open_box>;
    using periodic_searcher = md::neighbor_searcher_v2<md::periodic_box>;
    using xy_periodic_searcher = md::neighbor_searcher_v2<md::xy_periodic_box>;
    CHECK(sizeof(open_searcher) > 0);
    CHECK(sizeof(periodic_searcher) > 0);
    CHECK(sizeof(xy_periodic_searcher) > 0);
}
