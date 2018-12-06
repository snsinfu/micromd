#include <functional>
#include <memory>
#include <numeric>
#include <stdexcept>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/simulation/brownian_dynamics.hpp>

#include <catch.hpp>


TEST_CASE("brownian_dynamics_config - defaults to sane configuration")
{
    md::brownian_dynamics_config config;

    CHECK(config.temperature == 1);
    CHECK(config.timestep == 1);
    CHECK(config.steps == 1);
    CHECK(config.seed == 0);
    CHECK(!config.callback);
}

TEST_CASE("simulate_brownian_dynamics - does nothing if steps is zero")
{
    md::system system;

    system.add_particle();

    md::brownian_dynamics_config config;
    config.steps = 0;

    md::simulate_brownian_dynamics(system, config);

    CHECK(system.view_positions()[0].x == 0);
    CHECK(system.view_positions()[0].y == 0);
    CHECK(system.view_positions()[0].z == 0);
}

TEST_CASE("simulate_brownian_dynamics - can simulate frozen system")
{
    md::system system;

    system.add_particle();

    md::brownian_dynamics_config config;
    config.temperature = 0;
    config.steps = 100;

    md::simulate_brownian_dynamics(system, config);

    CHECK(system.view_positions()[0].x == 0);
    CHECK(system.view_positions()[0].y == 0);
    CHECK(system.view_positions()[0].z == 0);
}

TEST_CASE("simulate_brownian_dynamics - can simulate brownian motion")
{
    md::system system;

    system.add_particle();

    md::brownian_dynamics_config config;
    config.temperature = 1;
    config.steps = 100;

    md::simulate_brownian_dynamics(system, config);

    CHECK(system.view_positions()[0].x != 0);
    CHECK(system.view_positions()[0].y != 0);
    CHECK(system.view_positions()[0].z != 0);
}

TEST_CASE("simulate_brownian_dynamics - supports adaptive time-stepping")
{
    struct motion_recorder
    {
        md::system* system;
        md::scalar sum_motion = 0;
        md::step steps = 0;
        md::point prev_position;

        void operator()(md::step)
        {
            md::point const position = system->view_positions()[0];
            sum_motion += (position - prev_position).norm();
            prev_position = position;
            steps++;
        }

        md::scalar mean() const
        {
            return sum_motion / md::scalar(steps);
        }
    };

    md::system system;
    system.add_particle();
    system.add_particle();
    system.add_particle();

    motion_recorder recorder;
    recorder.system = &system;

    md::brownian_dynamics_config config;
    config.temperature = 1;
    config.spacestep = 0.01;
    config.steps = 1000;
    config.callback = std::ref(recorder);

    md::simulate_brownian_dynamics(system, config);

    CHECK(recorder.mean() == Approx(config.spacestep).epsilon(0.1));
}

TEST_CASE("simulate_brownian_dynamics - callback step is 1-based")
{
    md::system system;
    system.add_particle();

    std::vector<md::step> actual_steps;
    std::vector<md::step> expected_steps = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    md::brownian_dynamics_config config;
    config.steps = 10;
    config.callback = [&](md::step step) { actual_steps.push_back(step); };

    md::simulate_brownian_dynamics(system, config);

    CHECK(actual_steps == expected_steps);
}

TEST_CASE("simulate_brownian_dynamics - correctly samples simple canonical distribution")
{
    // Record radial position of a particle in a harmonic well.
    md::scalar const mu = 1.23;
    md::scalar const K = 4.56;
    md::scalar const kT = 7.89;
    md::scalar const dt = 0.98;
    md::step const sample_count = 5000;

    class test_forcefield : public md::forcefield
    {
    public:
        md::scalar K = 1;

        md::scalar compute_energy(md::system const&) override
        {
            throw std::logic_error("not implemented");
        }

        void compute_force(md::system const& system, md::array_view<md::vector> forces) override
        {
            md::array_view<md::point const> positions = system.view_positions();

            for (md::index i = 0; i < system.particle_count(); i++) {
                md::vector const r = positions[i] - md::point{0, 0, 0};
                md::vector const F = -K * r;
                forces[i] += F;
            }
        }
    };

    test_forcefield forcefield;
    forcefield.K = K;

    md::system system;
    system.add_particle().mobility = mu;
    system.add_forcefield(forcefield);

    // Sample radial positions from a Brownian trajectory.
    std::vector<md::scalar> radial_samples;

    auto const callback = [&](md::step) {
        md::point const pos = system.view_positions()[0];
        radial_samples.push_back(pos.distance(md::point{0, 0, 0}));
    };

    md::brownian_dynamics_config config;
    config.temperature = kT;
    config.timestep = dt;
    config.steps = sample_count;
    config.callback = callback;

    md::simulate_brownian_dynamics(system, config);

    // KS test (95%)
    md::scalar const critical_point = 1.63 / std::sqrt(sample_count);

    std::sort(radial_samples.begin(), radial_samples.end());

    std::vector<md::scalar> canonical_cdf;
    for (md::scalar const r : radial_samples) {
        md::scalar const pmf = r * r * std::exp(-0.5 * K * r * r / kT);
        canonical_cdf.push_back(pmf);
    }
    std::partial_sum(canonical_cdf.begin(), canonical_cdf.end(), canonical_cdf.begin());

    for (md::scalar& p : canonical_cdf) {
        p /= canonical_cdf.back();
    }

    md::scalar D = 0;
    md::scalar rank = 0;

    for (md::index i = 0; i < sample_count; i++) {
        md::scalar const sample_cdf = rank / md::scalar(sample_count);
        D = std::max(D, std::fabs(sample_cdf - canonical_cdf[i]));
        rank++;
    }

    CHECK(D < critical_point);
}
