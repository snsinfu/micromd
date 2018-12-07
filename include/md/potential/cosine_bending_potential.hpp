// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_POTENTIAL_COSINE_BENDING_POTENTIAL_HPP
#define MD_POTENTIAL_COSINE_BENDING_POTENTIAL_HPP

// This module provides the cosine bending potential.

#include <cmath>
#include <tuple>

#include "../basic_types.hpp"


namespace md
{
    // cosine_bending_potential implements the three-body bending potential
    // function:
    //
    //     u(r,s) = e dot(r,s) / |r||s|
    //
    struct cosine_bending_potential
    {
        md::scalar bending_energy = 0;

        md::scalar evaluate_energy(md::vector rij, md::vector rjk) const
        {
            md::scalar const dij_sq = rij.squared_norm();
            md::scalar const djk_sq = rjk.squared_norm();

            if (dij_sq * djk_sq == 0) {
                return 0;
            }

            md::scalar const dot_ijk = md::dot(rij, rjk);
            md::scalar const dij_djk = std::sqrt(dij_sq * djk_sq);
            md::scalar const cos_ijk = dot_ijk / dij_djk;

            return bending_energy * (1 - cos_ijk);
        }

        std::tuple<md::vector, md::vector, md::vector> evaluate_force(
            md::vector rij, md::vector rjk
        ) const
        {
            md::scalar const dij_sq = rij.squared_norm();
            md::scalar const djk_sq = rjk.squared_norm();

            if (dij_sq * djk_sq == 0) {
                return std::make_tuple(md::vector{}, md::vector{}, md::vector{});
            }

            md::scalar const dot_ijk = md::dot(rij, rjk);
            md::scalar const dij_djk = std::sqrt(dij_sq * djk_sq);
            md::scalar const e_div_dd = bending_energy / dij_djk;

            md::vector const fij = e_div_dd * (rjk - dot_ijk / dij_sq * rij);
            md::vector const fkj = e_div_dd * (rij - dot_ijk / djk_sq * rjk);

            return std::make_tuple(fij, -(fij + fkj), fkj);
        }
    };
}

#endif
