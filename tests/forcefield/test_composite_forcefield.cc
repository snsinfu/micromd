#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/forcefield/composite_forcefield.hpp>

#include <catch.hpp>


TEST_CASE("composite_forcefield - computes sum of composed forcefields")
{
    class my_component_a : public virtual md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const&) override
        {
            return 0.1;
        }

        void compute_force(md::system const&, md::array_view<md::vector> forces) override
        {
            for (md::vector& force : forces) {
                force.x += 0.1;
            }
        }
    };

    class my_component_b : public virtual md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const&) override
        {
            return 1;
        }

        void compute_force(md::system const&, md::array_view<md::vector> forces) override
        {
            for (md::vector& force : forces) {
                force.y += 1;
            }
        }
    };

    class my_forcefield : public md::composite_forcefield<my_component_a, my_component_b>
    {
    };

    md::system system;
    system.add_particle();

    my_forcefield forcefield;

    // Sum of component energies
    CHECK(forcefield.compute_energy(system) == Approx(1.1));

    // Sum of component forces
    std::vector<md::vector> forces(system.particle_count());
    forcefield.compute_force(system, forces);

    CHECK(forces[0].x == 0.1);
    CHECK(forces[0].y == 1);
    CHECK(forces[0].z == 0);
}

TEST_CASE("composite_forcefield - is zero forcefield when no component is given")
{
    class my_forcefield : public md::composite_forcefield<>
    {
    };

    md::system system;
    system.add_particle();

    my_forcefield forcefield;

    // Zero energy
    CHECK(forcefield.compute_energy(system) == 0);

    // Zero force
    std::vector<md::vector> forces(system.particle_count());
    forcefield.compute_force(system, forces);

    CHECK(forces[0].x == 0);
    CHECK(forces[0].y == 0);
    CHECK(forces[0].z == 0);
}

namespace
{
    template<typename Derived>
    class test_crtp_forcefield : public virtual md::forcefield
    {
    public:
        md::scalar compute_energy(md::system const&) override
        {
            return static_cast<Derived&>(*this).energy();
        }

        void compute_force(md::system const&, md::array_view<md::vector>) override
        {
        }
    };
}

TEST_CASE("composite_forcefield - is compatible with CRTP bases")
{
    class my_forcefield : public md::composite_forcefield<test_crtp_forcefield<my_forcefield>>
    {
    public:
        md::scalar energy() const
        {
            return 1.23;
        }
    };

    my_forcefield forcefield;
    md::system system;

    CHECK(forcefield.compute_energy(system) == 1.23);
}
