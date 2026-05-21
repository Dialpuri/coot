#pragma once
#include <clipper/core/coords.h>
#include <gemmi/model.hpp>

namespace coot {

void update_position_gemmi(gemmi::Atom *at, const clipper::Coord_orth &pos) {
   at->pos.x = pos.x();
   at->pos.y = pos.y();
   at->pos.z = pos.z();
}

} // namespace coot