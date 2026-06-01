#pragma once

#include <clipper/core/coords.h>
#include <gemmi/model.hpp>

namespace coot {
namespace util {

struct contact_atoms_info_t {
    struct contact_atom_t {
        double dist;
        gemmi::Atom const* at;
        double mat[4][4];

        // Factory that replaces the MMDB constructor
        static contact_atom_t contact_atom_t_gemmi(
            gemmi::Atom const* contactor,
            gemmi::Atom const* central_atom,
            double const mat_in[4][4])
        {
            contact_atom_t result;
            clipper::Coord_orth co_1(contactor->pos.x, contactor->pos.y, contactor->pos.z);
            clipper::Coord_orth co_2(central_atom->pos.x, central_atom->pos.y, central_atom->pos.z);
            result.dist = clipper::Coord_orth::length(co_1, co_2);
            result.at = contactor;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    result.mat[i][j] = mat_in[i][j];
            return result;
        }
    };
};

} // namespace util
} // namespace coot