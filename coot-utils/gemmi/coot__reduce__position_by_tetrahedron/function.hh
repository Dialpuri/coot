#pragma once

#include "clipper/core/coords.h"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"

namespace gemmi { class Atom; }

namespace coot { namespace reduce {

inline clipper::Coord_orth position_by_tetrahedron_gemmi(
    const gemmi::Atom* at_central,
    const gemmi::Atom* at_n_1,
    const gemmi::Atom* at_n_2,
    const gemmi::Atom* at_n_3,
    double bl) {

   clipper::Coord_orth at_c_pos = coot::co_gemmi(at_central);
   clipper::Coord_orth at_n1_pos = coot::co_gemmi(at_n_1);
   clipper::Coord_orth at_n2_pos = coot::co_gemmi(at_n_2);
   clipper::Coord_orth at_n3_pos = coot::co_gemmi(at_n_3);

   clipper::Coord_orth vec_1_to_c(at_c_pos - at_n1_pos);
   clipper::Coord_orth vec_2_to_c(at_c_pos - at_n2_pos);
   clipper::Coord_orth vec_3_to_c(at_c_pos - at_n3_pos);

   clipper::Coord_orth vec_1_to_c_uv(vec_1_to_c.unit());
   clipper::Coord_orth vec_2_to_c_uv(vec_2_to_c.unit());
   clipper::Coord_orth vec_3_to_c_uv(vec_3_to_c.unit());

   clipper::Coord_orth under_pos = 0.3333333333 *
      (vec_1_to_c_uv + vec_2_to_c_uv + vec_3_to_c_uv);
   clipper::Coord_orth under_pos_uv(under_pos.unit());

   clipper::Coord_orth H_pos = at_c_pos + bl * under_pos_uv;
   return H_pos;
}

}} // namespace coot::reduce