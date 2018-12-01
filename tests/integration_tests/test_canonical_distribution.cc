#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <vector>

#include <md/basic_types.hpp>
#include <md/system.hpp>
#include <md/forcefield/point_source_forcefield.hpp>
#include <md/potential/harmonic_potential.hpp>
#include <md/simulation/brownian_dynamics.hpp>

#include <catch.hpp>


TEST_CASE("Sampling canonical distribution with Brownian dynamics")
{
    // u(r) = Kr^2 sin(kr)
    // F(r) = -Kr ( 2 sin(kr) + kr cos(kr) )
    struct test_potential
    {
        md::scalar K = 1;
        md::scalar k = 2;

        md::scalar evaluate_energy(md::vector r) const
        {
            md::scalar const r1 = r.norm();
            md::scalar const r2 = r.squared_norm();

            return K * r2 * std::sin(k * r1);
        }

        md::vector evaluate_force(md::vector r) const
        {
            md::scalar const r1 = r.norm();
            md::scalar const kr = k * r1;

            return -K * (2 * std::sin(kr) + kr * std::cos(kr)) * r;
        }
    };

    class test_forcefield : public md::point_source_forcefield<test_forcefield>
    {
    public:
        test_potential point_source_potential(md::system const&, md::index) const
        {
            return test_potential{};
        }
    };

    // Sample radial positions
    md::system system;

    system.add_particle();
    system.add_forcefield(test_forcefield{});

    md::brownian_dynamics_config config;

    config.temperature = 1;
    config.timestep = 1e-5;
    config.spacestep = 0.01;
    config.steps = 1000000;

    std::vector<md::scalar> radial_samples;

    config.callback = [&](md::step step) {
        if ((step + 1) % 10 == 0) {
            md::point const pt = system.view_positions()[0];
            md::scalar const r = pt.vector().norm();
            radial_samples.push_back(r);
        }
    };

    md::simulate_brownian_dynamics(system, config);

    // Compute canonical distribution
    auto const canonical_density = [&](md::scalar r) {
        test_potential pot;
        md::scalar const u = pot.evaluate_energy(md::vector{r, 0, 0});
        return r * r * std::exp(-u / config.temperature);
    };

    std::vector<md::scalar> canonical_cdf;

    std::sort(radial_samples.begin(), radial_samples.end());

    std::transform(
        radial_samples.begin(),
        radial_samples.end(),
        std::back_inserter(canonical_cdf),
        canonical_density
    );

    std::partial_sum(canonical_cdf.begin(), canonical_cdf.end(), canonical_cdf.begin());

    for (md::scalar& dist : canonical_cdf) {
        dist /= canonical_cdf.back();
    }

    // KS test
    md::scalar const critical_value = 1.63 / std::sqrt(md::scalar(radial_samples.size()));

    md::scalar D = 0;

    for (std::size_t rank = 0; rank < radial_samples.size(); rank++) {
        md::scalar const F = md::scalar(rank) / md::scalar(radial_samples.size());
        D = std::max(D, std::fabs(F - canonical_cdf[rank]));
    }

    CHECK(D < critical_value);
}
