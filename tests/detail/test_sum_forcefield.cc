#include <memory>
#include <vector>

#include <md/basic_types.hpp>
#include <md/forcefield.hpp>
#include <md/system.hpp>

#include <md/detail/sum_forcefield.hpp>

#include "../catch.hpp"


TEST_CASE("sum_forcefield::add - adds a forcefield to the sum")
{
    class my_forcefield : public md::forcefield
    {
    public:
        md::scalar increase = 0;

        explicit my_forcefield(md::scalar increase)
            : increase{increase}
        {
        }

        md::scalar compute_energy(md::system const&) override
        {
            return increase;
        }

        void compute_force(md::system const&, md::array_view<md::vector> forces) override
        {
            for (md::vector& force : forces) {
                force.x += increase;
                force.y += increase;
                force.z += increase;
            }
        }
    };

    md::detail::sum_forcefield sum;

    sum.add(std::make_shared<my_forcefield>(1));
    sum.add(std::make_shared<my_forcefield>(2));

    md::system system;
    system.add_particle();

    // energy is the sum of component energy values
    CHECK(sum.compute_energy(system) == 1 + 2);

    // force is the sum of component force values
    std::vector<md::vector> forces(1);
    sum.compute_force(system, forces);
    CHECK(forces[0].x == 1 + 2);
    CHECK(forces[0].y == 1 + 2);
    CHECK(forces[0].z == 1 + 2);
}
