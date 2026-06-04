#pragma once

#include <gemmi/model.hpp>
#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace coot {

namespace reduce {

inline clipper::Coord_orth
position_by_bond_length_angle_torsion_gemmi(
    const gemmi::Atom *at_1,  // CA
    const gemmi::Atom *at_2,  // CB
    const gemmi::Atom *at_3,  // CG
    double bl,
    double angle_rad,
    double torsion_rad) {

   clipper::Coord_orth at_1_pos = coot::co_gemmi(at_1);
   clipper::Coord_orth at_2_pos = coot::co_gemmi(at_2);
   clipper::Coord_orth at_3_pos = coot::co_gemmi(at_3);

   clipper::Coord_orth new_pos(at_1_pos, at_2_pos, at_3_pos, bl, angle_rad, torsion_rad);

   return new_pos;
}

} // namespace reduce

} // namespace coot